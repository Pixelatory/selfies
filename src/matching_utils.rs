use std::{cmp::Reverse, collections::{BinaryHeap, HashMap, HashSet, VecDeque}};

use log::warn;


/// Finds a perfect matching for an undirected graph (without self-loops).
///
/// Graph: an adjacency list representing the input graph.
/// 
/// Returns a hashmap representing a perfect matching, where the key is the
/// graph node, and the value is the matching node. Returns None, if the graph
/// cannot be perfectly matched.
pub fn find_perfect_matching(graph: &HashMap<usize, Vec<usize>>) -> Option<HashMap<usize, usize>> {
    let mut matching = greedy_matching(graph)?;

    let mut unmatched: HashSet<usize> = (0..graph.len()).into_iter()
        .filter(|i| !matching.contains_key(i))
        .collect();

    while !unmatched.is_empty() {
        let first_item = *unmatched.iter().next()?;
        let root = unmatched.take(&first_item)?;
        
        let path = find_augmenting_path(graph, root, &matching)?;

        flip_augmenting_path(&mut matching, &path);
        unmatched.remove(path.first()?);
        unmatched.remove(path.last()?);
    }

    if matching.len() != graph.len() {
        // We should have already returned earlier if it was not
        // a perfect matching, but be defensive.
        None
    } else {
        Some(matching)
    }
}


/// Greedily match on the provided graph.
/// 
/// Picks the node with the least unmatched neighbours, and then matches
/// with the first unmatched neighbour. Repeats this process until no more
/// matchings can be made.
fn greedy_matching(graph: &HashMap<usize, Vec<usize>>) -> Option<HashMap<usize, usize>> {
    let mut matching = HashMap::new();
    // free_degrees[i] = number of unmatched neighbors for node i.
    let mut free_degrees: HashMap<usize, usize> = graph.iter().map(|(node, adj_nodes)| (*node, adj_nodes.len())).collect();

    // Prioritize nodes with fewer unmatched neighbours.
    let node_pqueue: Vec<Reverse<(usize, usize)>> = free_degrees.iter()
        .map(|(index, free_degree)| Reverse((*free_degree, *index)))
        .collect();
    let mut heap = BinaryHeap::from(node_pqueue);

    while !heap.is_empty() {
        let (_, node) = heap.pop()?.0;

        if matching.contains_key(&node) || *free_degrees.get(&node)? == 0 {
            // Node cannot be matched.
            continue;
        }
        
        // Match node with first unmatched neighbour.
        let mate = graph.get(&node)?
            .iter()
            .find(|x| !matching.contains_key(*x))
            .cloned()?;
        matching.insert(node, mate);
        matching.insert(mate, node);

        // Update the number of free degrees for neighbouring nodes.
        // Add them to the heap if they are not matched yet and have
        // candidate edges for matching.
        for adj in graph.get(&node)?.iter().chain(graph.get(&mate)?.iter()) {
            free_degrees.insert(*adj, free_degrees.get(adj)?.saturating_sub(1));
            if !matching.contains_key(adj) && *free_degrees.get(adj)? > 0 {
                heap.push(Reverse((free_degrees[adj], *adj)))
            }
        }
    }

    Some(matching)
}

/// Finds an augmenting path.
/// 
/// Root must be an unmatched node.
/// 
/// For matching problems, an augmenting path alternates between
/// traversing through a matched edge and an unmatched edge.
/// 
/// Returns the vertex path, or None if there is no augmenting path.
fn find_augmenting_path(graph: &HashMap<usize, Vec<usize>>, root: usize, matching: &HashMap<usize, usize>) -> Option<Vec<usize>> {
    if matching.contains_key(&root) {
        warn!("Attempted to find augmenting path with matched root.");
        return None;
    }

    // Run modified BFS to find path from root to unmatched node.
    let mut other_end = None;
    let mut node_queue = VecDeque::from(vec![root]);

    let mut parents: HashMap<usize, usize> = HashMap::new();

    while !node_queue.is_empty() {
        let node = node_queue.pop_front()?;

        for adj in graph.get(&node)?.iter().cloned() {
            if !matching.contains_key(&adj) && adj != root {
                parents.insert(adj, node);
                other_end = Some(adj);
                break;
            } else {
                let adj_mate =  *matching.get(&adj)?;
                if !parents.contains_key(&adj_mate) {
                    // adj_mate has not been visited.
                    parents.insert(adj_mate, adj);
                    parents.insert(adj, node);
                    node_queue.push_back(adj_mate);
                }
            }
        }

        if other_end.is_some() {
            // Augmenting path found.
            break;
        }
    }

    let mut path = vec![other_end?];
    let mut curr_node = other_end?;
    while curr_node != root {
        let parent = *parents.get(&curr_node)?;
        path.push(parent);
        curr_node = parent;
    }
    path.reverse();
    return Some(path);
}

/// Flips the matched nodes on the augmenting path.
/// 
/// ```text
/// I.e.
///          M
/// 0 --- 1 --- 2 --- 3
/// 
///    M           M
/// 0 --- 1 --- 2 --- 3
/// ```
fn flip_augmenting_path(matching: &mut HashMap<usize, usize>, path: &Vec<usize>) {
    for i in (0..path.len()).step_by(2) {
        let a = path.get(i);
        let b = path.get(i + 1);
        match (a, b) {
            (Some(aa), Some(bb)) => {
                matching.insert(*aa, *bb);
                matching.insert(*bb, *aa);
            },
            _ => continue,
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_augmenting_path() {
        // An augmenting path that traverses through matched nodes.
        let graph = HashMap::from(
            [
                (0, vec![1, 3]),
                (1, vec![0, 2]),
                (2, vec![1, 3]),
                (3, vec![2, 4]),
                (4, vec![3, 5]),
                (5, vec![4]),
            ]
        );
        let matching = HashMap::from([(1, 2), (2, 1), (3, 4), (4, 3)]);
        let x = find_augmenting_path(&graph, 0, &matching);
        assert!(x.is_some());
        assert_eq!(x.unwrap(), vec![0, 3, 4, 5]);

        // An augmenting path which goes from one unmatched node directly to another.
        let graph = HashMap::from(
            [
                (0, vec![1]),
                (1, vec![0, 2]),
                (2, vec![1, 3]),
                (3, vec![2, 4]),
                (4, vec![2, 3]),
            ]
        );
        let matching = HashMap::from([(1, 2),(2, 1)]);
        let x = find_augmenting_path(&graph, 4, &matching);
        assert!(x.is_some());
        assert_eq!(x.unwrap(), vec![4, 3]);
        
        // There is no augmenting path; no other unmatched node exists.
        let graph = HashMap::from(
            [
                (0, vec![1]),
                (1, vec![0, 2]),
                (2, vec![1]),
            ]
        );
        let matching = HashMap::from([(0, 1),(1, 0)]);
        let x = find_augmenting_path(&graph, 2, &matching);
        assert!(x.is_none());

        let graph = HashMap::from(
            [
                (0, vec![1]),
                (1, vec![2, 0, 3]),
                (2, vec![1]),
                (3, vec![1]),
            ]
        );
        let matching = HashMap::from([(0, 1),(1, 0)]);
        let x = find_augmenting_path(&graph, 3, &matching);
        assert!(x.is_none());
    }

    #[test]
    fn test_flip_augmenting_path() {
        // Flipping an empty matching.
        let mut matching = HashMap::new();
        let path = Vec::new();
        flip_augmenting_path(&mut matching, &path);
        assert_eq!(matching, HashMap::new());
        
        // Flip a valid matching.
        let mut matching = HashMap::from([(1, 2), (2, 1), (3, 4), (4, 3)]);
        let path = vec![0, 3, 4, 5];
        flip_augmenting_path(&mut matching, &path);
        // (1, 2), (2, 1) remain the same, however we flip out (3, 4), (4, 3).
        assert_eq!(matching, HashMap::from([(1, 2), (2, 1), (4, 5), (5, 4), (0, 3), (3, 0)]));
        
        // Flip when path goes from one unmatched node directly to another.
        let mut matching = HashMap::from([(1, 2),(2, 1)]);
        let path = vec![4, 3];
        flip_augmenting_path(&mut matching, &path);
        // (1, 2), (2, 1) remain, and we gain (3, 4), (4, 3).
        assert_eq!(matching, HashMap::from([(1, 2), (2, 1), (3, 4), (4, 3)]));
    }

    #[test]
    fn test_find_perfect_matching() {
        // A graph with a perfect matching.
        let graph = HashMap::from(
            [
                (0, vec![1, 3]),
                (1, vec![0, 2]),
                (2, vec![1, 3]),
                (3, vec![2, 4]),
                (4, vec![3, 5]),
                (5, vec![4]),
            ]
        );
        let x = find_perfect_matching(&graph);
        assert!(x.is_some());
        assert_eq!(x.unwrap(), HashMap::from([(0, 1), (1, 0), (2, 3), (3, 2), (4, 5), (5, 4)]));

        // A graph with no perfect matching.
        let graph = HashMap::from(
            [
                (0, vec![1]),
                (1, vec![2, 0, 3]),
                (2, vec![1]),
                (3, vec![1]),
            ]
        );
        let x = find_perfect_matching(&graph);
        assert!(x.is_none());
        
        // Another graph with no perfect matching.
        let graph = HashMap::from(
            [
                (0, vec![1]),
                (1, vec![0, 2]),
                (2, vec![1, 3]),
                (3, vec![2, 4]),
                (4, vec![2, 3]),
            ]
        );
        let x = find_perfect_matching(&graph);
        assert!(x.is_none());
        
        // A perfect matching on an empty graph is vacuously true.
        let graph = HashMap::new();
        let x = find_perfect_matching(&graph);
        assert!(x.is_some());
    }
}