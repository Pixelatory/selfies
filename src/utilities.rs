use crate::constants::ValenceTuple;


pub fn capitalize_first(s: &str) -> String {
    let mut c = s.chars();
    match c.next() {
        None => String::new(),
        Some(first) => first.to_uppercase().collect::<String>() + c.as_str(),
    }
}

pub fn last_valence(valence_tuple: &ValenceTuple) -> i32 {
    return match valence_tuple {
        ValenceTuple::One(x) => *x,
        ValenceTuple::Two(_, y) => *y,
    }
}

pub fn valence_any(valence_tuple: &ValenceTuple, check: i32) -> bool {
    return match valence_tuple {
        ValenceTuple::One(x) => *x == check,
        ValenceTuple::Two(x, y) => *x == check || *y == check,
    }
}