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

pub fn get_mut_pair<E>(collection: &mut [E], first_index: usize, second_index: usize) -> (&mut E, &mut E) {
    if second_index > first_index {
        let (first_half, second_half) = collection.split_at_mut(first_index + 1);
        return (&mut first_half[first_index], &mut second_half[second_index - first_index - 1]);
    } else {
        let (first_half, second_half) = collection.split_at_mut(second_index + 1);
        return (&mut second_half[first_index - second_index - 1], &mut first_half[second_index]);
    }
}