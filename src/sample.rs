use std::cmp;

extern crate rand;
use sample::rand::distributions::IndependentSample;

/// Sample k from n elements with replacement.
/// prob[i] represents the probability of sampling element i
pub fn sample(prob: &[f64], k: usize, mut rng: &mut rand::Rng) -> Vec<usize> {
    let n = prob.len();

    if n == 0 {
        return Vec::new();
    } else if n == 1 {
        return vec![0; k];
    }

    let mut p = indexed_revsorted(&prob);

    // convert probabilities into cumulative probabilities
    // which collectively represent a roulette
    for i in 1..n {
        p[i].0 += p[i-1].0;
    }

    // draw samples with replacement
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
        ys.push(p[j].1);
    }
   
    ys
}

/*
static void ProbSampleNoReplace(int n, double *p, int *perm,
				int nans, int *ans)
{
    double rT, mass, totalmass;
    int i, j, k, n1;

    /* Record element identities */
    for (i = 0; i < n; i++)
	perm[i] = i + 1;

    /* Sort probabilities into descending order */
    /* Order element identities in parallel */
    revsort(p, perm, n);

    /* Compute the sample */
    totalmass = 1;
    for (i = 0, n1 = n-1; i < nans; i++, n1--) {
	rT = totalmass * unif_rand();
	mass = 0;
	for (j = 0; j < n1; j++) {
	    mass += p[j];
	    if (rT <= mass)
		break;
	}
	ans[i] = perm[j];
	totalmass -= p[j];
	for(k = j; k < n1; k++) {
	    p[k] = p[k + 1];
	    perm[k] = perm[k + 1];
	}
    }
}
*/

/// Sample k from n elements without replacement.
/// prob[i] represents the probability of sampling element i
pub fn subsample(prob: &[f64], k: usize, mut rng: &mut rand::Rng) -> Vec<usize> {
    let n = prob.len();

    if n == 0 {
        return Vec::new();
    } else if n == 1 {
        return vec![0; k];
    } else if n <= k {
        return Vec::new();
    }

    let mut p = indexed_revsorted(&prob);

    // draw samples without replacment
    // i.e. previously drawn samples are removed
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

        let boot = sample(&prob, 10, &mut rng);
        assert_eq!(boot.len(), 10);

        let sub = subsample(&prob, 4, &mut rng);
        assert_eq!(sub.len(), 4);
    }
}
