pub fn flip_rows_in_place<T>(data: &mut [T], width: usize, height: usize) {
    if height <= 1 || width == 0 || data.is_empty() {
        return;
    }
    let row_len = width;
    for row in 0..(height / 2) {
        let top_start = row * row_len;
        let bottom_start = (height - 1 - row) * row_len;
        for col in 0..row_len {
            data.swap(top_start + col, bottom_start + col);
        }
    }
}

pub fn flip_rows_copy<T: Clone>(data: &[T], width: usize, height: usize) -> Vec<T> {
    if height <= 1 || width == 0 || data.is_empty() {
        return data.to_vec();
    }
    let mut result = Vec::with_capacity(data.len());
    for row in (0..height).rev() {
        let start = row * width;
        let end = start + width;
        result.extend_from_slice(&data[start..end]);
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flip_rows_in_place_2x2() {
        let mut data = vec![1, 2, 3, 4];
        flip_rows_in_place(&mut data, 2, 2);
        assert_eq!(data, vec![3, 4, 1, 2]);
    }

    #[test]
    fn flip_rows_in_place_3x3() {
        let mut data = vec![1, 2, 3, 4, 5, 6, 7, 8, 9];
        flip_rows_in_place(&mut data, 3, 3);
        assert_eq!(data, vec![7, 8, 9, 4, 5, 6, 1, 2, 3]);
    }

    #[test]
    fn flip_rows_in_place_1_row() {
        let mut data = vec![1, 2, 3];
        flip_rows_in_place(&mut data, 3, 1);
        assert_eq!(data, vec![1, 2, 3]);
    }

    #[test]
    fn flip_rows_in_place_empty() {
        let mut data: Vec<i32> = vec![];
        flip_rows_in_place(&mut data, 0, 0);
        assert!(data.is_empty());
    }

    #[test]
    fn flip_rows_copy_2x2() {
        let data = vec![1, 2, 3, 4];
        let result = flip_rows_copy(&data, 2, 2);
        assert_eq!(result, vec![3, 4, 1, 2]);
    }

    #[test]
    fn flip_rows_copy_preserves_original() {
        let data = vec![1, 2, 3, 4];
        let _result = flip_rows_copy(&data, 2, 2);
        assert_eq!(data, vec![1, 2, 3, 4]);
    }
}
