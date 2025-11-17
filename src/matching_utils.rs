
/// Finds a perfect matching for an undirected graph (without self-loops).
///
/// Graph: an adjacency list representing the input graph.
/// 
/// Returns a list representing a perfect matching, where j is the i-th
/// element if nodes i and j are matched. Returns None, if the graph cannot
/// be perfectly matched.
pub fn find_perfect_matching(graph: &Vec<Vec<usize>>) {
    let matching = greedy_matching(graph);
}

fn greedy_matching(graph: &Vec<Vec<usize>>) {

    // free_degrees[i] = number of unmatched neighbors for node i
    let free_degrees: Vec<usize> = graph.iter().map(|x| x.len()).collect();

    // prioritize nodes with fewer unmatched neighbours
    
}