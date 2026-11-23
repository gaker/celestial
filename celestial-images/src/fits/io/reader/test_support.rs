pub(super) struct MockFitsBuilder {
    cards: Vec<String>,
    data: Vec<u8>,
}

impl MockFitsBuilder {
    pub(super) fn new() -> Self {
        Self {
            cards: Vec::new(),
            data: Vec::new(),
        }
    }

    pub(super) fn card(mut self, keyword: &str, value: &str, comment: &str) -> Self {
        let card = if comment.is_empty() {
            format!("{:<8}= {:<70}", keyword, value)
        } else {
            format!("{:<8}= {:<20} / {}", keyword, value, comment)
        };
        let mut card_80 = format!("{:<80}", card);
        card_80.truncate(80);
        self.cards.push(card_80);
        self
    }

    fn simple_primary(self) -> Self {
        self.card("SIMPLE", "T", "Standard FITS format")
            .card("BITPIX", "8", "Bits per pixel")
            .card("NAXIS", "0", "Number of axes")
            .card("EXTEND", "F", "No extensions")
    }

    fn image_with_data<T>(self, bitpix: i32, dims: &[usize], data: &[T]) -> Self
    where
        T: Copy + Into<f64>,
    {
        let mut builder = self
            .card("SIMPLE", "T", "Standard FITS format")
            .card("BITPIX", &bitpix.to_string(), "Bits per pixel")
            .card("NAXIS", &dims.len().to_string(), "Number of axes");

        for (i, &dim) in dims.iter().enumerate() {
            let keyword = format!("NAXIS{}", i + 1);
            builder = builder.card(&keyword, &dim.to_string(), &format!("Axis {} size", i + 1));
        }

        builder = builder.card("EXTEND", "F", "No extensions");

        let mut builder_with_data = MockFitsBuilder {
            cards: builder.cards,
            data: Vec::new(),
        };

        match bitpix {
            8 => {
                for &val in data {
                    builder_with_data.data.push(val.into() as u8);
                }
            }
            16 => {
                for &val in data {
                    let bytes = (val.into() as i16).to_be_bytes();
                    builder_with_data.data.extend_from_slice(&bytes);
                }
            }
            32 => {
                for &val in data {
                    let bytes = (val.into() as i32).to_be_bytes();
                    builder_with_data.data.extend_from_slice(&bytes);
                }
            }
            -32 => {
                for &val in data {
                    let bytes = (val.into() as f32).to_be_bytes();
                    builder_with_data.data.extend_from_slice(&bytes);
                }
            }
            -64 => {
                for &val in data {
                    let bytes = val.into().to_be_bytes();
                    builder_with_data.data.extend_from_slice(&bytes);
                }
            }
            _ => panic!("Unsupported BITPIX: {}", bitpix),
        }

        builder_with_data
    }

    pub(super) fn build_memory(mut self) -> Vec<u8> {
        self.cards.push("END".to_string() + &" ".repeat(77));

        let mut fits_data = Vec::new();

        for card in &self.cards {
            let mut card_bytes = card.as_bytes().to_vec();
            card_bytes.resize(80, b' ');
            fits_data.extend_from_slice(&card_bytes);
        }

        while fits_data.len() % 2880 != 0 {
            fits_data.push(b' ');
        }

        if !self.data.is_empty() {
            fits_data.extend_from_slice(&self.data);

            while fits_data.len() % 2880 != 0 {
                fits_data.push(0);
            }
        }

        fits_data
    }
}

pub(super) fn create_minimal_fits() -> Vec<u8> {
    MockFitsBuilder::new().simple_primary().build_memory()
}

pub(super) fn create_image_fits<T>(bitpix: i32, dims: &[usize], data: &[T]) -> Vec<u8>
where
    T: Copy + Into<f64>,
{
    MockFitsBuilder::new()
        .image_with_data(bitpix, dims, data)
        .build_memory()
}

pub(super) fn create_malformed_fits(issue: &str) -> Vec<u8> {
    match issue {
        "no_end" => {
            let mut data = Vec::new();
            let cards = [
                format!(
                    "{:<8}= {:<20} / {:<47}",
                    "SIMPLE", "T", "Standard FITS format"
                ),
                format!("{:<8}= {:<20} / {:<47}", "BITPIX", "8", "Bits per pixel"),
                format!("{:<8}= {:<20} / {:<47}", "NAXIS", "0", "Number of axes"),
            ];

            for card in &cards {
                let mut card_bytes = card.as_bytes().to_vec();
                card_bytes.resize(80, b' ');
                data.extend_from_slice(&card_bytes);
            }
            data.resize(2880, b' ');
            data
        }
        "invalid_utf8" => {
            let mut data = vec![0xFF, 0xFE, 0xFD];
            data.resize(2880, b' ');
            data
        }
        "truncated" => {
            vec![0x53, 0x49, 0x4D]
        }
        _ => create_minimal_fits(),
    }
}
