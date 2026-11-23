// Axis lengths joined by "x", as in 7x5x3, for error messages.
pub(crate) fn shape_text(axes: &[usize]) -> String {
    let lengths: Vec<String> = axes.iter().map(usize::to_string).collect();
    lengths.join("x")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shape_text_joins_axis_lengths() {
        assert_eq!(shape_text(&[7, 5, 3]), "7x5x3");
        assert_eq!(shape_text(&[9]), "9");
    }
}
