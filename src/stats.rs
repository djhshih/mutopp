use std::f64;
use std::cmp;

/// Calculate log(sum(exp(xs))) of a real vector xs.
/// NaN values are ignored.
pub fn log_sum_exp(xs: &[f64]) -> f64 {
    // find maximum value
    let mut x_max = f64::NAN;
    for &x in xs.iter() {
        if x > x_max || x_max.is_nan() {
            x_max = x;
        }
    }

    if x_max.is_nan() {
        return f64::NAN;
    } else if x_max.is_infinite() {
        return x_max;
    }

    // sum the differences
    let mut sum = 0.0;
    for &x in xs.iter() {
        if !x.is_nan() {
            sum += (x - x_max).exp();
        }
    }

    x_max + sum.ln()
}

/// Vector dot product
pub fn dot_product(xs: &[f64], ys: &[f64]) -> f64 {
    let n = cmp::min(xs.len(), ys.len());
    let xs = &xs[..n];
    let ys = &ys[..n];

    let mut s = 0.0;
    for i in 0..n {
        s += xs[i] * ys[i];
    }

    s
}

/// Cosine similarity
/// bounded in [-1, 1]
pub fn cosine_similarity(xs: &[f64], ys: &[f64]) -> f64 {
   dot_product(xs, ys) / (dot_product(xs, xs) * dot_product(ys, ys)).sqrt()
}

pub fn scale(xs: &[f64]) -> Vec<f64> {
    let s: f64 = xs.iter().sum();
    if s == 1.0 {
        xs.to_owned()
    } else {
        xs.iter().map(|x| x / s).collect()
    }
}

/// Kullback-Leibler Divergence using log2
pub fn kl_divergence(p: &[f64], q: &[f64]) -> f64 {
    let n = cmp::min(p.len(), q.len());
    let p = &p[..n];
    let q = &q[..n];

    let mut d = 0.0;
    for i in 0..n {
        if p[i] != 0.0 {
            d += p[i] * (p[i].log2() - q[i].log2())
        }
    }
    
    d
}

/// Jensen-Shannon Divergence using log2
/// bounded in [0, 1]
pub fn js_divergence(p: &[f64], q: &[f64]) -> f64 {
    // ensure that p and q are proper distributions
    let p = scale(p);
    let q = scale(q);

    let m: Vec<f64> = p.iter().zip(q.iter()).map(|(&x, &y)| 0.5 * (x + y)).collect();
    
    0.5 * (kl_divergence(&p, &m) + kl_divergence(&q, &m))
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1.0e-6;

    #[test]
    fn test_log_sum_exp() {
        assert!((log_sum_exp(&vec![0.2, 0.3, 0.9, 1.2]) - 2.122379).abs() < EPS);
        assert!((log_sum_exp(&vec![0.0, 1.0, -5.0]) - 1.315072).abs() < EPS);
        assert!((log_sum_exp(&vec![0.0, 1.0, f64::NAN]) - 1.313262).abs() < EPS);
        assert!(log_sum_exp(&vec![0.0, f64::NEG_INFINITY]) == 0.0);
        assert!(log_sum_exp(&vec![0.0, f64::INFINITY]).is_infinite());
    }

    #[test]
    fn test_js_divergence() {
        assert!(js_divergence(&vec![0.0, 0.1, 0.9], &vec![0.0, 0.1, 0.9]) == 0.0);
        assert!((js_divergence(&vec![0.0, 0.1, 0.9], &vec![0.1, 0.2, 0.7]) - 0.0712961).abs() < EPS);
        assert!((js_divergence(&vec![0.1, 0.2, 0.7], &vec![0.0, 0.1, 0.9]) - 0.0712961).abs() < EPS);
    }

    #[test]
    fn test_kl_divergence() {
        assert!(kl_divergence(&vec![0.0, 0.1, 0.9], &vec![0.0, 0.1, 0.9]) == 0.0);
        assert!((kl_divergence(&vec![0.1, 0.1, 0.8], &vec![0.1, 0.2, 0.7]) - 0.05411606).abs() < EPS);
        assert!((kl_divergence(&vec![0.1, 0.2, 0.7], &vec![0.1, 0.1, 0.8]) - 0.06514845).abs() < EPS);
    }

    #[test]
    fn test_cosine_similarity() {
        assert!((cosine_similarity(&vec![0.0, 0.1, 0.9], &vec![0.0, 0.1, 0.9]) - 1.0).abs() < EPS);
        assert!((cosine_similarity(&vec![0.1, 0.1, 0.8], &vec![0.1, 0.2, 0.7]) - 0.9882872).abs() < EPS);
        assert!((cosine_similarity(&vec![0.1, 0.2, 0.7], &vec![0.1, 0.1, 0.8]) - 0.9882872).abs() < EPS);
    }
}
