use std::cmp::{Ordering, PartialOrd};

/// A point within a file.
/// Zero-based, so the first row is zero.
#[derive(Debug, Default, PartialEq, Eq, PartialOrd, Ord, Copy, Clone, Hash)]
pub struct Point {
    pub row: usize,
    pub column: usize,
}

impl From<tree_sitter::Point> for Point {
    fn from(value: tree_sitter::Point) -> Self {
        Point {
            row: value.row,
            column: value.column,
        }
    }
}

impl From<lsp_types::Position> for Point {
    fn from(value: lsp_types::Position) -> Self {
        Point {
            row: value.line as usize,
            column: value.character as usize,
        }
    }
}

impl Point {
    /// Converts into a `tree_sitter::Point`
    pub fn into_tree_sitter(self) -> tree_sitter::Point {
        tree_sitter::Point {
            row: self.row,
            column: self.column,
        }
    }

    /// Converts into an `lsp_types::Position`
    /// Note: may cause issues if row/column > `u32::MAX`
    pub fn into_lsp_type(self) -> lsp_types::Position {
        lsp_types::Position {
            line: self.row as u32,
            character: self.column as u32,
        }
    }
}

impl From<(usize, usize)> for Point {
    fn from((row, column): (usize, usize)) -> Self {
        Self { row, column }
    }
}

#[derive(Debug, Default, PartialEq, Eq, PartialOrd, Ord, Copy, Clone, Hash)]
pub struct Span {
    pub start: Point,
    pub end: Point,
}

impl From<tree_sitter::Range> for Span {
    fn from(value: tree_sitter::Range) -> Self {
        Self {
            start: value.start_point.into(),
            end: value.end_point.into(),
        }
    }
}

impl From<lsp_types::Range> for Span {
    fn from(value: lsp_types::Range) -> Self {
        Self {
            start: value.start.into(),
            end: value.end.into(),
        }
    }
}

impl From<((usize, usize), (usize, usize))> for Span {
    fn from((start, end): ((usize, usize), (usize, usize))) -> Self {
        let start = start.into();
        let end = end.into();

        Self { start, end }
    }
}
impl From<(Point, Point)> for Span {
    fn from((start, end): (Point, Point)) -> Self {
        Self { start, end }
    }
}

impl<P> From<std::ops::Range<P>> for Span
where
    P: Into<Point>,
{
    fn from(value: std::ops::Range<P>) -> Self {
        Self {
            start: value.start.into(),
            end: value.end.into(),
        }
    }
}

impl Span {
    /// Convert the span into a `tree-sitter::Range`.
    /// Requires additional information about the bytes span of the token to be provided.
    pub fn into_tree_sitter(self, start_byte: usize, end_byte: usize) -> tree_sitter::Range {
        tree_sitter::Range {
            start_point: self.start.into_tree_sitter(),
            end_point: self.end.into_tree_sitter(),
            start_byte,
            end_byte,
        }
    }

    /// Convert the span into an `lsp_types::Range`
    /// NOTE: if the line/column number is > `u32::MAX`, casting from `usize` may result in issues.
    pub fn into_lsp_types(self) -> lsp_types::Range {
        lsp_types::Range {
            start: self.start.into_lsp_type(),
            end: self.end.into_lsp_type(),
        }
    }
}

impl Span {
    pub fn in_range(&self, range: &Span) -> bool {
        range.start <= self.start && self.end <= range.end
    }

    /// Whether the point is in the spans range, inclusively
    pub fn contains(&self, point: &Point) -> bool {
        self.start <= *point && *point <= self.end
    }

    /// Whether the point is in the spans range, inclusively
    pub fn contains_span(&self, other: &Span) -> bool {
        self.start <= other.start && other.end <= self.end
    }
}

impl PartialOrd<Point> for Span {
    /// Compare a point and a span
    ///
    /// If point is inside the span considered equal.
    /// If outside to the left, the span is considered greater than the point
    /// If outside to the right, considered greater.
    /// Both ends of the range are inclusive
    fn partial_cmp(&self, other: &Point) -> Option<std::cmp::Ordering> {
        match (self.start.cmp(other), self.end.cmp(other)) {
            // End is less than the point, span must be less than point
            (_, Ordering::Less) => Some(Ordering::Less),
            // Start is greater than the poinbt, span must be greater than point
            (Ordering::Greater, _) => Some(Ordering::Greater),
            // Means the point is inside the span.
            (Ordering::Less, Ordering::Greater) => Some(Ordering::Equal),
            (_, Ordering::Equal) => Some(Ordering::Equal),
            (Ordering::Equal, _) => Some(Ordering::Equal),
        }
    }
}

impl PartialOrd<Span> for Point {
    /// Compare a point and a span
    ///
    /// If point is inside the span considered equal.
    /// If outside to the left, the span is considered greater than the point
    /// If outside to the right, considered greater.
    /// Both ends of the range are inclusive
    fn partial_cmp(&self, other: &Span) -> Option<std::cmp::Ordering> {
        match (self.cmp(&other.start), self.cmp(&other.end)) {
            // Point is less than the start, must be outside
            (Ordering::Less, _) => Some(Ordering::Less),
            // Point is greater than the start, must be outside
            (_, Ordering::Greater) => Some(Ordering::Greater),
            // TODO: Explicitly enumerate
            (Ordering::Greater, Ordering::Less) => Some(Ordering::Equal),
            (Ordering::Equal, _) => Some(Ordering::Equal),
            (_, Ordering::Equal) => Some(Ordering::Equal),
        }
    }
}

impl PartialEq<Point> for Span {
    /// Compare a point and a span
    ///
    /// If point is inside the span considered equal.
    /// If outside to the left, considered less than.
    /// If outside to the right, considered greater.
    /// Both ends of the range are inclusive
    fn eq(&self, other: &Point) -> bool {
        self.contains(other)
    }
}

impl PartialEq<Span> for Point {
    /// Compare a point and a span
    ///
    /// If point is inside the span considered equal.
    /// If outside to the left, considered less than.
    /// If outside to the right, considered greater.
    /// Both ends of the range are inclusive
    fn eq(&self, other: &Span) -> bool {
        other.contains(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn span_contains() {
        let point = Point { row: 0, column: 2 };
        let span = Span {
            start: Point { row: 0, column: 1 },
            end: Point { row: 0, column: 3 },
        };

        assert!(span.contains(&point));
        assert_eq!(span, point);

        assert_eq!(span.partial_cmp(&point), Some(Ordering::Equal));
    }

    #[test]
    fn point_on_lower_edge_of_span() {
        let point = Point { row: 0, column: 1 };
        let span = Span {
            start: Point { row: 0, column: 1 },
            end: Point { row: 0, column: 3 },
        };

        assert!(span.contains(&point));
        assert_eq!(span, point);

        assert_eq!(span.partial_cmp(&point), Some(Ordering::Equal));
    }

    #[test]
    fn point_on_outside_left_of_span() {
        let point = Point { row: 0, column: 0 };
        let span = Span {
            start: Point { row: 0, column: 1 },
            end: Point { row: 0, column: 3 },
        };

        assert!(!span.contains(&point));
        assert!(point < span);
        assert!(span > point);

        assert_eq!(point.partial_cmp(&span), Some(Ordering::Less));
        assert_eq!(span.partial_cmp(&point), Some(Ordering::Greater));
    }

    #[test]
    fn point_on_outside_right_of_span() {
        let point = Point { row: 0, column: 10 };
        let span = Span {
            start: Point { row: 0, column: 1 },
            end: Point { row: 0, column: 3 },
        };

        assert!(!span.contains(&point));
        assert!(point > span);
        assert!(span < point);

        assert_eq!(point.partial_cmp(&span), Some(Ordering::Greater));
        assert_eq!(span.partial_cmp(&point), Some(Ordering::Less));
    }

    #[test]
    fn point_on_upper_edge_of_span() {
        let point = Point { row: 0, column: 3 };
        let span = Span {
            start: Point { row: 0, column: 1 },
            end: Point { row: 0, column: 3 },
        };

        assert!(span.contains(&point));
        assert_eq!(span, point);

        assert_eq!(span.partial_cmp(&point), Some(Ordering::Equal));
    }

    #[test]
    fn point_ordering() {
        let p0 = Point { row: 0, column: 0 };
        let p1 = Point { row: 1, column: 0 };

        assert!(p1 > p0);

        let p2 = Point { row: 0, column: 1 };
        let p3 = Point { row: 0, column: 2 };

        assert!(p2 > p0);
        assert!(p3 > p0);

        // Higher row always higher order
        assert!(p1 > p2);
        assert!(p1 > p3);
    }
}
