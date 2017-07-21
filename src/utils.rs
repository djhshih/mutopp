/// Parse genomic coordinate.
///
/// Assume coordinate is specified in 1-based, closed format.
/// Internally, coordinate is converted to 0-based, half-open format.
///
pub fn parse_coordinate(x: &str) -> Option<(String, u64, u64)> {
    match x.find(':') {
        Some(i) => {
            let seqname = &x[..i];
            let range = &x[i+1..];
            match range.find('-') {
                Some(j) => {
                    let start = &range[..j];
                    let end = &range[j+1..];
                    if let Ok(start) = start.parse() {
                        if let Ok(end) = end.parse() {
                            if start > 0 && end > start {
                                Some((String::from(seqname), start-1, end))
                            } else {
                                None
                            }
                        } else {
                            None
                        }
                    } else {
                        None
                    }
                },
                None => None,
            }
        },
        None => None,
    }
}