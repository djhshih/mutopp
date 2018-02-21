extern crate rand;

/// Sample k from n elements without replacement.
pub fn subsample(n: usize, k: usize, prob: &[u64], rng: &mut rand::Rng) -> Vec<usize> {

    // make indexed probability tuple
    let mut p: Vec<(u64, usize)> = Vec::new();
    for i in 0..prob.len() {
        p[i] = (prob[i], i);
    }
   
    return Vec::new();
}

/// Sort vector xs by descending order using heapsort,
/// and sort vector idx by the same permuation.
/// If idx == 0..n initially, idx will contain the sorted permutation.
// idx must be the same size as xs
fn revsort<T: Copy + PartialOrd>(xs: &mut [T], idx: &mut [usize]) {
    let n = xs.len();

    if n <= 1 {
        return;
    }
    
    let mut l = n / 2;
    let mut ir = n - 1;

    let mut ra: T;
    let mut ii: usize;

    let mut i;
    let mut j;

    loop {
        if l > 0 {
            l = l - 1;
            ra = xs[l];
            ii = idx[l];
        } else {
            ra = xs[ir];
            ii = idx[ir];

            xs[ir] = xs[0];
            idx[ir] = idx[0];

            ir -= 1;
            if ir == 0 {
                xs[0] = ra;
                idx[0] = ii;
                return;
            }
        }

        i = l;
        j = l / 2;
        while j <= ir {
            if j < ir && xs[j] > xs[j + 1] {
                j += 1;
            }
            if ra > xs[j] {
                xs[i] = xs[j];
                idx[i] = idx[j];
                i = j;
                j += j;
            } else {
                j = ir + 1;
            }
        }
        xs[i] = ra;
        idx[i] = ii;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_revsort() {
        let mut xs = vec![1.0, 5.2, 3.4, 3.3];
        let mut idx = vec![0, 1, 2, 3];

        revsort(&mut xs, &mut idx);

        println!("{:?}", xs);
        println!("{:?}", idx);

        assert_eq!(xs, vec![1.0, 3.3, 3.4, 5.2]);
    }
}
