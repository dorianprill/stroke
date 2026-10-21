//! Subdivision model and shared rotation search.
//!
//! The parent module documents the method, references and error limits.

use super::*;

const MAX_LEAVES: usize = 64;
const INITIAL_LEAVES: usize = 16;
const MAX_SEEDS: usize = 8;
const ACTIVE_SEEDS: usize = 4;
const MODEL_TRIALS_PER_SEED: usize = 512;
// Target the model gap at selected frames. This does not bound optimization error.
const MODEL_RELATIVE_VOLUME_TARGET: f64 = 1e-3;

// Cache subcurve controls so further splits need no extraction from the original.
// Endpoint projections lower-bound widths. Control-hull projections upper-bound
// them in exact arithmetic. Floating-point rounding is not certified.
struct Model<T, const D: usize, const N: usize, const L: usize = MAX_LEAVES> {
    pieces: [[[T; D]; N]; L],
    // Flattening nested arrays gives room for both endpoints at jumps without
    // imposing another generic-const-expression bound on the model.
    points: [[[T; D]; 2]; L],
    point_count: usize,
    parameters: [[T; 2]; L],
    errors: [T; L],
    blocked: [bool; L],
    len: usize,
}

fn split<T: Float, const D: usize, const N: usize>(
    mut points: [[T; D]; N],
    t: T,
) -> ([[T; D]; N], [[T; D]; N]) {
    let mut left = points;
    let mut right = points;
    let s = T::one() - t;
    for level in 1..N {
        for i in 0..N - level {
            points[i] = array::from_fn(|d| s * points[i][d] + t * points[i + 1][d]);
        }
        left[level] = points[0];
        right[N - level - 1] = points[N - level - 1];
    }
    (left, right)
}

fn chord_error<T: Float, const D: usize, const N: usize>(points: &[[T; D]; N]) -> T {
    let first = points[0];
    let chord: [T; D] = array::from_fn(|d| points[N - 1][d] - first[d]);
    let length2 = dot(&chord, &chord);
    points.iter().fold(T::zero(), |error, point| {
        let offset = array::from_fn(|d| point[d] - first[d]);
        let t = if length2 > T::zero() {
            (dot(&offset, &chord) / length2)
                .max(T::zero())
                .min(T::one())
        } else {
            T::zero()
        };
        let v: [T; D] = array::from_fn(|d| (-t).mul_add(chord[d], offset[d]));
        error.max(dot(&v, &v).sqrt())
    })
}

fn sample_bounds<T: Float, const D: usize>(points: &[[T; D]], axes: &Frame<T, D>) -> [(T, T); D] {
    axes.map(|axis| {
        points
            .iter()
            .fold((T::infinity(), T::neg_infinity()), |(lo, hi), p| {
                let x = dot(&axis, p);
                (lo.min(x), hi.max(x))
            })
    })
}

impl<T: Float, const D: usize, const N: usize, const L: usize> Model<T, D, N, L> {
    fn new(controls: [[T; D]; N]) -> Self {
        let mut model = Self {
            pieces: [controls; L],
            points: [[[T::zero(); D]; 2]; L],
            point_count: 0,
            parameters: [[T::zero(), T::one()]; L],
            errors: [T::zero(); L],
            blocked: [false; L],
            len: 1,
        };
        model.rebuild_points();
        model.errors[0] = chord_error(&controls);
        model
    }

    fn from_pieces(
        mut pieces: impl Iterator<Item = [[T; D]; N]>,
    ) -> Result<Self, OrientedBoundingBoxError> {
        let first = pieces
            .next()
            .ok_or(OrientedBoundingBoxError::NonFiniteComputation)?;
        let mut model = Self::new(first);
        for piece in pieces {
            if model.len == L {
                return Err(OrientedBoundingBoxError::NonFiniteComputation);
            }
            model.pieces[model.len] = piece;
            model.errors[model.len] = chord_error(&piece);
            model.len += 1;
        }
        if model.pieces[..model.len]
            .iter()
            .flatten()
            .flatten()
            .any(|x| !x.is_finite())
        {
            return Err(OrientedBoundingBoxError::NonFiniteComputation);
        }
        model.rebuild_points();
        model.refine_flatness(INITIAL_LEAVES.min(L));
        Ok(model)
    }

    fn rebuild_points(&mut self) {
        let points = self.points.as_flattened_mut();
        self.point_count = 0;
        for piece in &self.pieces[..self.len] {
            for point in [piece[0], piece[N - 1]] {
                if self.point_count == 0 || points[self.point_count - 1] != point {
                    points[self.point_count] = point;
                    self.point_count += 1;
                }
            }
        }
    }

    fn points(&self) -> &[[T; D]] {
        &self.points.as_flattened()[..self.point_count]
    }

    fn error(&self) -> T {
        self.errors[..self.len]
            .iter()
            .copied()
            .fold(T::zero(), T::max)
    }

    fn subdivide(&mut self, index: usize) -> bool {
        if self.len == L || self.blocked[index] {
            return false;
        }
        let [a, b] = self.parameters[index];
        let middle = a * scalar::<T>(0.5) + b * scalar::<T>(0.5);
        if middle <= a || middle >= b {
            self.blocked[index] = true;
            return false;
        }
        let (left, right) = split(self.pieces[index], scalar(0.5));
        for i in (index + 1..self.len).rev() {
            self.pieces[i + 1] = self.pieces[i];
            self.errors[i + 1] = self.errors[i];
            self.blocked[i + 1] = self.blocked[i];
            self.parameters[i + 1] = self.parameters[i];
        }
        self.pieces[index] = left;
        self.pieces[index + 1] = right;
        self.errors[index] = chord_error(&left);
        self.errors[index + 1] = chord_error(&right);
        self.blocked[index] = false;
        self.blocked[index + 1] = false;
        self.parameters[index] = [a, middle];
        self.parameters[index + 1] = [middle, b];
        self.len += 1;
        self.rebuild_points();
        true
    }

    fn refine_flatness(&mut self, limit: usize) {
        let tolerance = scalar::<T>(1e-4).max(scalar::<T>(128.0) * T::epsilon());
        while self.len < limit {
            let index = (0..self.len)
                .filter(|&i| !self.blocked[i])
                .max_by(|&a, &b| {
                    self.errors[a]
                        .partial_cmp(&self.errors[b])
                        .unwrap_or(core::cmp::Ordering::Equal)
                });
            let Some(index) = index else { break };
            if self.errors[index] <= tolerance {
                break;
            }
            self.subdivide(index);
        }
    }

    fn upper_bounds(&self, axes: &Frame<T, D>) -> [(T, T); D] {
        axes.map(|axis| {
            self.pieces[..self.len].iter().flatten().fold(
                (T::infinity(), T::neg_infinity()),
                |(lo, hi), p| {
                    let x = dot(&axis, p);
                    (lo.min(x), hi.max(x))
                },
            )
        })
    }

    // Split the piece whose controls most extend an uncertain box face, relative
    // to that side's width. Pieces inside the endpoint bounds need no split.
    fn refine_frames<const C: usize>(
        &mut self,
        frames: &[Option<Candidate<T, D>>; MAX_SEEDS],
        data: &Normalized<T, D, C>,
    ) {
        let tolerance = scalar::<T>(MODEL_RELATIVE_VOLUME_TARGET).ln_1p();
        while self.len < L {
            let mut best = None;
            let mut largest = T::zero();
            for candidate in frames.iter().flatten() {
                if self.log_gap(&candidate.axes, data) <= tolerance {
                    continue;
                }
                let lower = sample_bounds(self.points(), &candidate.axes);
                for (d, axis) in candidate.axes.iter().enumerate() {
                    let padding = frame_padding(data, axis);
                    let width = lower[d].1 - lower[d].0 + scalar::<T>(2.0) * padding;
                    for i in 0..self.len {
                        if self.blocked[i] {
                            continue;
                        }
                        let (lo, hi) = self.pieces[i].iter().fold(
                            (T::infinity(), T::neg_infinity()),
                            |(lo, hi), p| {
                                let x = dot(axis, p);
                                (lo.min(x), hi.max(x))
                            },
                        );
                        let excess = (lower[d].0 - lo).max(hi - lower[d].1).max(T::zero());
                        let relative = excess / width;
                        if relative > largest {
                            largest = relative;
                            best = Some(i);
                        }
                    }
                }
            }
            let Some(index) = best else { break };
            self.subdivide(index);
        }
    }

    fn log_gap<const C: usize>(&self, axes: &Frame<T, D>, data: &Normalized<T, D, C>) -> T {
        let lower = sample_bounds(self.points(), axes);
        let upper = self.upper_bounds(axes);
        (0..D).fold(T::zero(), |sum, i| {
            let padding = scalar::<T>(2.0) * frame_padding(data, &axes[i]);
            let lo = lower[i].1 - lower[i].0 + padding;
            let hi = upper[i].1 - upper[i].0 + padding;
            sum + (hi / lo).max(T::one()).ln()
        })
    }
}

fn frame_padding<T: Float, const D: usize, const N: usize>(
    data: &Normalized<T, D, N>,
    axis: &[T; D],
) -> T {
    axis.iter()
        .zip(data.roundoff)
        .fold(data.margin, |sum, (&a, e)| a.abs().mul_add(e, sum))
}

fn residual<T: Float, const D: usize>(mut v: [T; D], axes: &[[T; D]]) -> [T; D] {
    for _ in 0..2 {
        for axis in axes {
            let p = dot(&v, axis);
            for (x, &a) in v.iter_mut().zip(axis) {
                *x = (-p).mul_add(a, *x);
            }
        }
    }
    v
}

fn geometric_frame<T: Float, const D: usize>(
    points: &[[T; D]],
    a: usize,
    b: usize,
) -> Option<Frame<T, D>> {
    let mut axes = identity();
    axes[0] = normalize(array::from_fn(|d| points[b][d] - points[a][d]))?;
    let floor = scalar::<T>(128.0) * T::epsilon() * <T as NumCast>::from(D).unwrap();
    for i in 1..D {
        let mut best = [T::zero(); D];
        let mut length2 = T::zero();
        for point in points {
            let v = residual(array::from_fn(|d| point[d] - points[a][d]), &axes[..i]);
            let n = dot(&v, &v);
            if n > length2 {
                best = v;
                length2 = n;
            }
        }
        if length2 <= floor * floor {
            length2 = T::zero();
            for d in 0..D {
                let v = residual(
                    array::from_fn(|j| if d == j { T::one() } else { T::zero() }),
                    &axes[..i],
                );
                let n = dot(&v, &v);
                if n > length2 {
                    best = v;
                    length2 = n;
                }
            }
        }
        axes[i] = normalize(best)?;
    }
    Some(axes)
}

struct Seeds<T, const D: usize> {
    entries: [Option<Candidate<T, D>>; MAX_SEEDS],
    diameter: T,
}
impl<T: Float, const D: usize> Seeds<T, D> {
    fn new(diameter: T) -> Self {
        Self {
            entries: [None; MAX_SEEDS],
            diameter,
        }
    }
    fn similar(&self, a: &Candidate<T, D>, b: &Candidate<T, D>) -> bool {
        let wa = a.bounds.extent();
        let wb = b.bounds.extent();
        a.axes.iter().enumerate().all(|(i, axis)| {
            b.axes.iter().enumerate().any(|(j, other)| {
                let sign = T::one().copysign(dot(axis, other));
                let v: [T; D] = array::from_fn(|d| axis[d] - sign * other[d]);
                let tolerance = (scalar::<T>(0.05) * wa[i].min(wb[j]) / self.diameter)
                    .max(scalar::<T>(32.0) * T::epsilon())
                    .min(scalar(0.1));
                dot(&v, &v) <= tolerance * tolerance
            })
        })
    }
    fn insert(&mut self, candidate: Candidate<T, D>) {
        if let Some(i) = self
            .entries
            .iter()
            .position(|e| e.is_some_and(|other| self.similar(&candidate, &other)))
        {
            if !improves(&candidate, &self.entries[i].unwrap()) {
                return;
            }
            self.entries[i..].rotate_left(1);
            self.entries[MAX_SEEDS - 1] = None;
        }
        if let Some(i) = self
            .entries
            .iter()
            .position(|e| e.is_none_or(|other| improves(&candidate, &other)))
        {
            for j in (i + 1..MAX_SEEDS).rev() {
                self.entries[j] = self.entries[j - 1];
            }
            self.entries[i] = Some(candidate);
        }
    }
}

// A poll tests small rotations around the current frame. Keep its cursor
// between trials so every budget follows the same sequence. Combined rotations
// can improve a box when separate rotations stall at a change of extreme point.
struct Poll<T, const D: usize> {
    best: Candidate<T, D>,
    base: Candidate<T, D>,
    step: T,
    cursor: usize,
    cycle: u64,
    changed: bool,
    done: bool,
}
impl<T: Float, const D: usize> Poll<T, D> {
    fn new(best: Candidate<T, D>, step: T) -> Self {
        Self {
            best,
            base: best,
            step,
            cursor: 0,
            cycle: 0,
            changed: false,
            done: false,
        }
    }
    fn next<const N: usize>(&mut self, data: &Normalized<T, D, N>) -> Option<Frame<T, D>> {
        let pairs = D * (D - 1) / 2;
        let coordinate_trials = 2 * pairs;
        let total = coordinate_trials + if pairs > 1 { 4 } else { 0 };
        if self.cursor == total {
            if !self.changed {
                self.step = self.step * scalar(0.5);
            }
            self.cursor = 0;
            self.cycle = self.cycle.wrapping_add(1);
        }
        if self.cursor == 0 {
            self.base = self.best;
            self.changed = false;
            let widths = self.base.bounds.extent();
            let noise = (0..D).fold(T::zero(), |sum, i| {
                sum + scalar::<T>(2.0) * frame_padding(data, &self.base.axes[i]) / widths[i]
            });
            let floor = (noise * scalar(0.25))
                .max(scalar::<T>(16.0) * T::epsilon())
                .min(scalar(0.25));
            if self.step < floor {
                self.done = true;
                return None;
            }
        }
        let widths = self.base.bounds.extent();
        let mut frame = self.base.axes;
        let mut index = 0;
        for i in 0..D {
            for j in i + 1..D {
                let coefficient = if self.cursor < coordinate_trials {
                    if index == self.cursor / 2 {
                        if self.cursor % 2 == 0 {
                            -T::one()
                        } else {
                            T::one()
                        }
                    } else {
                        T::zero()
                    }
                } else {
                    let joint = self.cursor - coordinate_trials;
                    let sign = if joint % 2 == 0 { -T::one() } else { T::one() };
                    let value = if joint < 2 {
                        T::one()
                    } else {
                        let mut hash = self
                            .cycle
                            .wrapping_add(1)
                            .wrapping_mul(0x9e37_79b9_7f4a_7c15)
                            .wrapping_add(index as u64);
                        hash = (hash ^ (hash >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
                        hash = (hash ^ (hash >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
                        hash ^= hash >> 31;
                        scalar::<T>((hash >> 11) as f64 / ((1_u64 << 53) as f64) * 2.0 - 1.0)
                    };
                    sign * value / <T as NumCast>::from(pairs).unwrap().sqrt()
                };
                if coefficient != T::zero() {
                    let ratio = widths[i].min(widths[j]) / widths[i].max(widths[j]);
                    frame = rotated(frame, i, j, self.step * ratio * coefficient);
                }
                index += 1;
            }
        }
        self.cursor += 1;
        Some(frame)
    }
    fn accept(&mut self, trial: Candidate<T, D>) {
        if improves(&trial, &self.best) {
            self.best = trial;
            self.changed = true;
        }
    }
}

pub(crate) fn fit<T: Float, P: PointIndex<Scalar = T>, const N: usize>(
    curve: &Bezier<P, N>,
    budget: usize,
) -> Result<OrientedBoundingBox<T, { P::DIM }>, OrientedBoundingBoxError> {
    check_request(P::DIM, budget)?;
    let original =
        array::from_fn::<_, N, _>(|i| array::from_fn(|d| curve.control_points().as_slice()[i][d]));
    let data = normalized_controls(
        original,
        curve.control_box().center(),
        <T as NumCast>::from(N).unwrap(),
    )?;
    if budget == 1 || P::DIM == 1 || data.scale == T::zero() {
        return initial_box(&data, &mut |axes| projected_bounds(&data.points, axes));
    }
    search(
        &data,
        budget,
        Model::<T, { P::DIM }, N>::from_pieces(core::iter::once(data.points))?,
        &mut |axes| projected_bounds(&data.points, axes),
    )
}

pub(crate) fn fit_bspline<
    T: Float,
    P: PointIndex<Scalar = T>,
    const K: usize,
    const C: usize,
    const D: usize,
>(
    curve: &BSpline<P, K, C, D>,
    budget: usize,
) -> Result<OrientedBoundingBox<T, { P::DIM }>, OrientedBoundingBoxError>
where
    [(); D + 1]: Sized,
    [(); C + MAX_LEAVES]: Sized,
{
    check_request(P::DIM, budget)?;
    let (data, normalized) = normalized_bspline(curve)?;
    if budget == 1 || P::DIM == 1 || data.scale == T::zero() {
        return initial_box(&data, &mut |axes| {
            projected_bspline_bounds(&normalized, axes)
        });
    }
    let model = Model::<T, { P::DIM }, { D + 1 }, { C + MAX_LEAVES }>::from_pieces(
        normalized
            .bezier_spans()
            .map(|(_, span)| span.into_control_points().map(ArrayPoint::into_array)),
    )?;
    search(&data, budget, model, &mut |axes| {
        projected_bspline_bounds(&normalized, axes)
    })
}

fn search<T: Float, const D: usize, const C: usize, const N: usize, const L: usize>(
    data: &Normalized<T, D, C>,
    budget: usize,
    mut model: Model<T, D, N, L>,
    bounds_for_frame: &mut impl FnMut(&Frame<T, D>) -> Result<[(T, T); D], OrientedBoundingBoxError>,
) -> Result<OrientedBoundingBox<T, D>, OrientedBoundingBoxError> {
    // Callers handle constants, 1D and budget one before building the model.
    let mut best = evaluate(&data.roundoff, data.margin, identity(), bounds_for_frame)?;
    let mut evaluations = 1;
    let diameter = data
        .points
        .iter()
        .fold([T::zero(); D], |mut widths, p| {
            for d in 0..D {
                widths[d] = widths[d].max((p[d] - data.points[0][d]).abs());
            }
            widths
        })
        .iter()
        .fold(T::zero(), |sum, &w| sum + w * w)
        .sqrt()
        * scalar::<T>(2.0);
    let mut seeds = Seeds::new(diameter);
    let mut model_evaluations = 0;
    let mut evaluate_model = |frame| {
        model_evaluations += 1;
        evaluate(&data.roundoff, data.margin, frame, &mut |axes| {
            Ok(model.upper_bounds(axes))
        })
    };
    seeds.insert(evaluate_model(identity())?);
    // Limit starting chords for splines with many spans. Complete and score
    // each frame using all model points and pieces.
    let count = model.points().len().min(INITIAL_LEAVES + 1);
    for i in 0..count {
        for j in i + 1..count {
            let a = i * (model.points().len() - 1) / (count - 1);
            let b = j * (model.points().len() - 1) / (count - 1);
            if let Some(frame) = geometric_frame(model.points(), a, b) {
                seeds.insert(evaluate_model(frame)?);
            }
        }
    }
    // Reserve independent starts so similar frames do not fill the whole queue.
    // Nearby orientations can lead to different minima when extreme points change.
    let first = model.points()[0];
    let chord = array::from_fn(|d| model.points()[model.points().len() - 1][d] - first[d]);
    let farthest = model
        .points()
        .iter()
        .map(|p| array::from_fn(|d| p[d] - first[d]))
        .max_by(|a, b| {
            dot(a, a)
                .partial_cmp(&dot(b, b))
                .unwrap_or(core::cmp::Ordering::Equal)
        })
        .unwrap();
    for (index, frame) in [
        Some(identity()),
        aligned_frame(chord),
        aligned_frame(farthest),
    ]
    .into_iter()
    .enumerate()
    {
        if let Some(frame) = frame {
            seeds.entries[MAX_SEEDS - 3 + index] = Some(evaluate_model(frame)?);
        }
    }
    model.refine_frames(&seeds.entries, data);
    let mut states: [Option<Poll<T, D>>; MAX_SEEDS] = array::from_fn(|_| None);
    for (i, seed) in seeds.entries.into_iter().flatten().enumerate() {
        let mut state = Poll::new(
            evaluate(&data.roundoff, data.margin, seed.axes, &mut |axes| {
                Ok(model.upper_bounds(axes))
            })?,
            scalar(0.5),
        );
        model_evaluations += 1;
        for _ in 0..MODEL_TRIALS_PER_SEED {
            let Some(frame) = state.next(data) else {
                break;
            };
            let trial = evaluate(&data.roundoff, data.margin, frame, &mut |axes| {
                Ok(model.upper_bounds(axes))
            })?;
            model_evaluations += 1;
            state.accept(trial);
        }
        states[i] = Some(state);
    }
    let frames = array::from_fn(|i| states[i].as_ref().map(|s| s.best));
    model.refine_frames(&frames, data);
    let mut refined = Seeds::new(diameter);
    for state in states.iter().flatten() {
        refined.insert(evaluate(
            &data.roundoff,
            data.margin,
            state.best.axes,
            &mut |axes| Ok(model.upper_bounds(axes)),
        )?);
        model_evaluations += 1;
    }
    let mut active: [Option<Poll<T, D>>; ACTIVE_SEEDS] = array::from_fn(|_| None);
    for (i, seed) in refined
        .entries
        .into_iter()
        .flatten()
        .take(ACTIVE_SEEDS)
        .enumerate()
    {
        if evaluations == budget {
            break;
        }
        let trial = evaluate(&data.roundoff, data.margin, seed.axes, bounds_for_frame)?;
        evaluations += 1;
        if improves(&trial, &best) {
            best = trial;
        }
        let previous = states
            .iter()
            .flatten()
            .find(|s| refined.similar(&s.best, &seed))
            .map_or(scalar(0.01), |s| s.step);
        let step = (previous * scalar(2.0))
            .max(model.log_gap(&seed.axes, data))
            .min(scalar(0.5));
        active[i] = Some(Poll::new(trial, step));
    }
    while evaluations < budget {
        let mut worked = false;
        for state in active.iter_mut().flatten() {
            if evaluations == budget {
                break;
            }
            if state.done {
                continue;
            }
            let Some(frame) = state.next(data) else {
                continue;
            };
            let trial = evaluate(&data.roundoff, data.margin, frame, bounds_for_frame)?;
            evaluations += 1;
            worked = true;
            state.accept(trial);
            if improves(&trial, &best) {
                best = trial;
            }
        }
        if !worked {
            break;
        }
    }
    let diagnostics = OrientedBoundingBoxDiagnostics {
        representative_points: model.points().len(),
        representative_limit_reached: model.len == L
            || model.blocked[..model.len].iter().any(|&b| b),
        point_set_evaluations: model_evaluations,
        model_log_volume_gap: model.log_gap(&best.axes, data),
        representative_error: model.error() * data.scale,
    };
    Ok(OrientedBoundingBox {
        origin: data.origin,
        scale: data.scale,
        candidate: best,
        evaluations,
        diagnostics: Some(diagnostics),
    })
}
#[cfg(test)]
mod tests {
    use super::super::tests::cubic_range;
    use super::*;

    fn representatives<T: Float, const D: usize, const N: usize>(
        controls: [[T; D]; N],
    ) -> Model<T, D, N> {
        let mut model = Model::new(controls);
        model.refine_flatness(MAX_LEAVES);
        model
    }

    fn random(state: &mut u64) -> f64 {
        *state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*state >> 11) as f64 / ((1_u64 << 53) as f64)) * 2.0 - 1.0
    }

    fn check_widths<const D: usize>(controls: [[f64; D]; 4], state: &mut u64) {
        let samples = representatives(controls);
        assert!((1..=MAX_LEAVES).contains(&samples.len));
        assert!(samples.error().is_finite() && samples.error() >= 0.0);
        assert_eq!(samples.points()[0], controls[0]);
        assert_eq!(samples.points().last(), Some(&controls[3]));
        for _ in 0..64 {
            let axis = normalize(array::from_fn::<_, D, _>(|_| random(state))).unwrap();
            let (lo, hi) = cubic_range(controls.map(|p| dot(&axis, &p)));
            let bounds = sample_bounds(samples.points(), &array::from_fn(|_| axis));
            let (a, b) = bounds[0];
            let roundoff = 256.0 * D as f64 * f64::EPSILON;
            // These checks use analytic extrema, not a parameter grid.
            assert!(lo >= a - samples.error() - roundoff);
            assert!(hi <= b + samples.error() + roundoff);
            assert!(a >= lo - roundoff && b <= hi + roundoff);
            assert!((hi - lo) - (b - a) <= 2.0 * samples.error() + 2.0 * roundoff);
        }
    }

    #[test]
    fn additive_directional_width_bound_in_2d_3d_4d() {
        let mut state = 0x000c_05e7;
        for _ in 0..32 {
            check_widths(
                array::from_fn::<_, 4, _>(|_| array::from_fn::<_, 2, _>(|_| random(&mut state))),
                &mut state,
            );
            check_widths(
                array::from_fn::<_, 4, _>(|_| array::from_fn::<_, 3, _>(|_| random(&mut state))),
                &mut state,
            );
            check_widths(
                array::from_fn::<_, 4, _>(|_| array::from_fn::<_, 4, _>(|_| random(&mut state))),
                &mut state,
            );
        }
    }

    #[test]
    fn disconnected_piece_model_brackets_all_directional_widths() {
        let pieces = [
            [
                [-1.0_f64, 0.0, 0.1],
                [-0.8, 1.0, 0.5],
                [-0.5, -0.8, 0.0],
                [-0.3, 0.0, 0.3],
            ],
            [
                [0.6, 0.4, -0.2],
                [0.8, -0.6, 0.5],
                [0.7, 0.5, 0.3],
                [1.0, -0.2, 0.1],
            ],
        ];
        let mut model = Model::<_, 3, 4, 80>::from_pieces(pieces.into_iter()).unwrap();
        model.refine_flatness(80);
        assert!(model.points().contains(&pieces[0][3]));
        assert!(model.points().contains(&pieces[1][0]));
        let mut state = 0xb51;
        for _ in 0..64 {
            let axis = normalize(array::from_fn(|_| random(&mut state))).unwrap();
            let (lo, hi) = pieces
                .iter()
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| {
                    let (a, b) = cubic_range(p.map(|p| dot(&axis, &p)));
                    (lo.min(a), hi.max(b))
                });
            let frame = [axis; 3];
            let (a, b) = sample_bounds(model.points(), &frame)[0];
            let (outer_lo, outer_hi) = model.upper_bounds(&frame)[0];
            assert!(a >= lo - 1e-12 && b <= hi + 1e-12);
            assert!(outer_lo <= lo + 1e-12 && outer_hi >= hi - 1e-12);
            assert!((hi - lo) - (b - a) <= 2.0 * model.error() + 1e-12);
        }
    }

    #[test]
    fn representative_capacity_retains_error_and_handles_backtracking() {
        let controls = [[-1.0_f64, -1.0], [1.0, 1.0], [-1.0, 1.0], [1.0, -1.0]];
        let samples = representatives(controls);
        assert_eq!(samples.len, MAX_LEAVES);
        assert!(samples.error() > 1e-4);
        let mut state = 0x51;
        check_widths(controls, &mut state);
        check_widths([[0.0_f64], [4.0], [-3.0], [1.0]], &mut state);
        check_widths(
            [[0.0_f64; 3], [1.0, 1.0, 1.0], [-1.0, 1.0, -1.0], [0.0; 3]],
            &mut state,
        );
        let line = representatives([[0.0_f64, 0.0], [1.0, 2.0], [2.0, 4.0], [3.0, 6.0]]);
        assert_eq!(line.len, 1);
        assert!(line.error() < 32.0 * f64::EPSILON);
    }

    #[test]
    fn projected_hulls_bracket_extrema_and_tighten_on_subdivision() {
        let controls = [
            [0.0_f64, 0.0, 0.0],
            [1.0, 3.0, 2.0],
            [4.0, -2.0, 1.0],
            [5.0, 1.0, 4.0],
        ];
        let mut model = Model::<_, 3, 4>::from_pieces(core::iter::once(controls)).unwrap();
        let mut state = 0xabba;
        let axes =
            orthogonalize(array::from_fn(|_| array::from_fn(|_| random(&mut state)))).unwrap();
        let lower = sample_bounds(model.points(), &axes);
        let upper = model.upper_bounds(&axes);
        model.subdivide(0);
        let refined_lower = sample_bounds(model.points(), &axes);
        let refined_upper = model.upper_bounds(&axes);
        for i in 0..3 {
            let (lo, hi) = cubic_range(controls.map(|p| dot(&axes[i], &p)));
            assert!(upper[i].0 <= lo + 1e-12 && upper[i].1 >= hi - 1e-12);
            assert!(lower[i].0 >= lo - 1e-12 && lower[i].1 <= hi + 1e-12);
            assert!(
                refined_upper[i].0 >= upper[i].0 - 1e-12
                    && refined_upper[i].1 <= upper[i].1 + 1e-12
            );
            assert!(
                refined_lower[i].0 <= lower[i].0 + 1e-12
                    && refined_lower[i].1 >= lower[i].1 - 1e-12
            );
        }
        let data = normalized_controls(controls, [0.0; 3], 4.0).unwrap();
        let candidate = evaluate(&data.roundoff, data.margin, axes, &mut |a| {
            Ok(model.upper_bounds(a))
        })
        .unwrap();
        let gap = model.log_gap(&axes, &data);
        let mut frames = [None; MAX_SEEDS];
        frames[0] = Some(candidate);
        model.refine_frames(&frames, &data);
        assert!(model.log_gap(&axes, &data) <= gap + 1e-12);
        if model.len < MAX_LEAVES && !model.blocked[..model.len].iter().any(|&b| b) {
            assert!(model.log_gap(&axes, &data) <= MODEL_RELATIVE_VOLUME_TARGET.ln_1p() + 1e-12);
        }
    }

    #[test]
    fn parameter_resolution_stops_subdivision_without_losing_the_hull() {
        let mut model = Model::<_, 1, 2>::new([[0.0_f64], [1.0]]);
        model.parameters[0][0] = 1.0 - f64::EPSILON / 2.0;
        let before = model.upper_bounds(&identity());
        assert!(!model.subdivide(0));
        assert!(model.blocked[0]);
        assert_eq!(model.len, 1);
        assert_eq!(model.upper_bounds(&identity()), before);
    }

    #[test]
    fn coordinated_poll_escapes_a_coordinatewise_stall() {
        let data = normalized_controls([[0.0_f64; 3], [1.0; 3]], [0.0; 3], 2.0).unwrap();
        let initial = Candidate {
            axes: identity(),
            bounds: BoundingBox::new([-1.0; 3], [1.0; 3]).unwrap(),
            score: (0.0, 0.0),
        };
        let mut poll = Poll::new(initial, 0.1);
        for _ in 0..6 {
            let axes = poll.next(&data).unwrap();
            let objective = axes[0][1].max(axes[0][2]);
            assert!(objective >= 0.0);
            poll.accept(Candidate {
                axes,
                score: (objective, 0.0),
                ..initial
            });
        }
        let axes = poll.next(&data).unwrap();
        let objective = axes[0][1].max(axes[0][2]);
        assert!(objective < 0.0);
        poll.accept(Candidate {
            axes,
            score: (objective, 0.0),
            ..initial
        });
        assert!(poll.best.score.0 < 0.0);
    }

    #[test]
    fn thin_frames_are_not_deduplicated_by_a_fixed_angle() {
        let bounds = BoundingBox::new([0.0; 3], [1.0, 0.1, 1e-6]).unwrap();
        let a = Candidate {
            axes: identity(),
            bounds,
            score: (0.0, 0.0),
        };
        let b = Candidate {
            axes: rotated(identity(), 0, 2, 1e-5),
            ..a
        };
        let seeds = Seeds::new(2.0);
        assert!(!seeds.similar(&a, &b));
        let same = Candidate {
            axes: [a.axes[2].map(|x| -x), a.axes[0], a.axes[1]],
            bounds: BoundingBox::new([0.0; 3], [1e-6, 1.0, 0.1]).unwrap(),
            ..a
        };
        assert!(seeds.similar(&a, &same));
    }
}
