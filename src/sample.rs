extern crate rand;
use sample::rand::distributions::IndependentSample;

/// Sample k from n elements.
/// prob[i] represents the probability of sampling element i
/// prob[i] >= 0 && sum(prob) = 1
pub fn sample(prob: &[f64], k: usize, replace: bool, mut rng: &mut rand::Rng) -> Vec<usize> {
    if replace {
        sample_replace(prob, k, &mut rng)
    } else {
        sample_no_replace(prob, k, &mut rng)
    }
}

/// Resample elements s.t. the distribution of a factor variable matches a target distribution,
/// using rejection sampling.
/// factor[i] contains factor value for observation i
/// target is the target distribution for the factor
/// min_freq is the minimum observed frequency to match
pub fn resample_factor_distrib(factor: &[usize], target: &[f64], min_freq: f64, mut rng: &mut rand::Rng) -> Vec<usize> {
    // number of observations
    let n =  factor.len();
    // number of factor levels
    let p = target.len();
    
    // calculate the observed distribution of factor values
    let mut observed = vec![0 as f64; p];
    for &f in factor.iter() {
        observed[f] += 1.0;
    }
    let total: f64 = observed.iter().sum();
    for x in observed.iter_mut() {
        *x /= total;
    }

    // calculate value m s.t. target[j] / (m * observed[j]) < 1 for all j
    // the acceptance probability is 1/m
    let mut m = 0.0;
    for j in 0..p {
        // if any factor frequency has very low frequency, there is little
        // we can do about it, since we cannot draw new samples;
        // therefore, we do the best we can to match frequencies of other factor values
        if observed[j] > min_freq {
            let s = target[j] / observed[j];
            if s > m {
                m = s;
            }
        }
    }

    let mut ys: Vec<usize> = Vec::with_capacity(n);
    if m > 0.0 {
        let r = rand::distributions::Range::new(0.0, 1.0);
        for (i, &f) in factor.iter().enumerate() {
            let u = r.ind_sample(&mut rng);
            if u < target[f] / (observed[f] * m) {
                // accept sample
                ys.push(i);
            }
        }
    }

    ys
}

/// Sample k from n elements with replacement.
// TODO replace with more efficient algorithm
fn sample_replace(prob: &[f64], k: usize, mut rng: &mut rand::Rng) -> Vec<usize> {
    let n = prob.len();

    if n == 0 {
        return Vec::new();
    } else if n == 1 {
        return vec![0; k];
    }

    // time O(n log n)
    let mut p = indexed_revsorted(&prob);

    // convert probabilities into cumulative probabilities
    // which collectively represent a roulette
    for i in 1..n {
        p[i].0 += p[i-1].0;
    }

    // draw samples with replacement
    // time O(n k)
    // i.e. previously drawn samples are not removed
    let mut ys = Vec::with_capacity(k);
    let r = rand::distributions::Range::new(0.0, 1.0);
    let last = n - 1;
    for _ in 0..k {
        let u = r.ind_sample(&mut rng);
        // find index j where random value u lands
        let mut j = 0;
        while j < last && u > p[j].0 {
            j += 1;
        }
        // record the element index
        ys.push(p[j].1);
    }
   
    ys
}

/// Sample k from n elements without replacement.
// TODO consider sampling with replacement but either
// 1. keep track of sampled elements using a map, or
// 2. remove duplicates later
fn sample_no_replace(prob: &[f64], k: usize, mut rng: &mut rand::Rng) -> Vec<usize> {
    let n = prob.len();

    if n == 0 || n < k {
        return Vec::new();
    } else if n == 1 {
        return vec![0; k];
    }

    // time O(n log n)
    let mut p = indexed_revsorted(&prob);

    // draw samples without replacment
    // i.e. previously drawn samples are removed
    // time O(k n)
    let mut ys = Vec::with_capacity(k);
    let r = rand::distributions::Range::new(0.0, 1.0);
    let mut total = 1.0;
    let mut last = n - 1;
    for _ in 0..k {
        let t = total * r.ind_sample(&mut rng);

        // find index j where random value t lands;
        // mass needs to be dynamically computed
        // because p is mutable
        let mut mass = p[0].0;
        let mut j = 0;
        while j < last && t > mass {
            mass += p[j].0;
            j += 1;
        }
        // record the element index
        ys.push(p[j].1);

        // remove element j from the roulette
        // by shifting elements left
        for k in j..last {
            p[k] = p[k+1];
        }
        total -= p[j].0;
        last -= 1;
    }
   
    ys
}

/// Sort vector xs by descending order with its index.
fn indexed_revsorted(xs: &[f64]) -> Vec<(f64, usize)>  {
    if xs.is_empty() {
        return Vec::new();
    }

    // make (probability, index) tuples
    let mut ts: Vec<(f64, usize)> = Vec::with_capacity(xs.len());
    for (i, x) in xs.iter().enumerate() {
        ts.push((*x, i));
    }

    // sort the proabilities in descending order
    ts.sort_unstable_by(|a, b| b.0.partial_cmp(&a.0).unwrap());

    ts
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sort() {
        let xs = vec![1.0, 5.2, 3.4, 3.3];
        let ys = indexed_revsorted(&xs);

        let idx: Vec<usize> = ys.iter().map(|&x| x.1).collect();

        assert_eq!(idx, vec![1, 2, 3, 0]);
    }

    #[test]
    fn test_sample() {
        let prob = vec![0.5, 0.2, 0.15, 0.1, 0.05];
        
        let mut rng = rand::thread_rng();

        let boot = sample(&prob, 10, true, &mut rng);
        assert_eq!(boot.len(), 10);

        let sub = sample(&prob, 4, false, &mut rng);
        assert_eq!(sub.len(), 4);

        let invalid = sample(&prob, 20, false, &mut rng);
        assert!(invalid.is_empty());
    }
}
