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

/*
/* Unequal probability sampling; with-replacement case */

static void ProbSampleReplace(int n, double *p, int *perm, int nans, int *ans)
{
    double rU;
    int i, j;
    int nm1 = n - 1;

    /* record element identities */
    for (i = 0; i < n; i++)
	perm[i] = i + 1;

    /* sort the probabilities into descending order */
    revsort(p, perm, n);

    /* compute cumulative probabilities */
    for (i = 1 ; i < n; i++)
	p[i] += p[i - 1];

    /* compute the sample */
    for (i = 0; i < nans; i++) {
	rU = unif_rand();
	for (j = 0; j < nm1; j++) {
	    if (rU <= p[j])
		break;
	}
	ans[i] = perm[j];
    }
}
*/

/*
void revsort(double *a, int *ib, int n)
{
/* Sort a[] into descending order by "heapsort";
 * sort ib[] alongside;
 * if initially, ib[] = 1...n, it will contain the permutation finally
 */

    int l, j, ir, i;
    double ra;
    int ii;

    if (n <= 1) return;

    a--; ib--;

    l = (n >> 1) + 1;
    ir = n;

    for (;;) {
	if (l > 1) {
	    l = l - 1;
	    ra = a[l];
	    ii = ib[l];
	}
	else {
	    ra = a[ir];
	    ii = ib[ir];
	    a[ir] = a[1];
	    ib[ir] = ib[1];
	    if (--ir == 1) {
		a[1] = ra;
		ib[1] = ii;
		return;
	    }
	}
	i = l;
	j = l << 1;
	while (j <= ir) {
	    if (j < ir && a[j] > a[j + 1]) ++j;
	    if (ra > a[j]) {
		a[i] = a[j];
		ib[i] = ib[j];
		j += (i = j);
	    }
	    else
		j = ir + 1;
	}
	a[i] = ra;
	ib[i] = ii;
    }
}
*/

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
