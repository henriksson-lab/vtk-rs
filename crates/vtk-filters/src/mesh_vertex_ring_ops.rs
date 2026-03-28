//! Vertex ring operations: 1-ring, n-ring extraction and analysis.

use vtk_data::{AnyDataArray, DataArray, Points, PolyData};

/// Extract the 1-ring neighborhood of a vertex.
pub fn one_ring(mesh: &PolyData, vertex: usize) -> Vec<usize> {
    let n=mesh.points.len();
    let adj=build_adj(mesh,n);
    if vertex<n{adj[vertex].clone()}else{Vec::new()}
}

/// Extract the n-ring neighborhood (all vertices within n hops).
pub fn n_ring(mesh: &PolyData, vertex: usize, radius: usize) -> Vec<usize> {
    let n=mesh.points.len();
    let adj=build_adj(mesh,n);
    let mut visited=std::collections::HashSet::new();
    let mut queue=std::collections::VecDeque::new();
    if vertex<n{queue.push_back((vertex,0));visited.insert(vertex);}
    while let Some((v,d))=queue.pop_front(){
        if d>=radius{continue;}
        for &nb in &adj[v]{if visited.insert(nb){queue.push_back((nb,d+1));}}
    }
    visited.into_iter().collect()
}

/// Compute per-vertex valence (number of adjacent vertices).
pub fn vertex_valence(mesh: &PolyData) -> PolyData {
    let n=mesh.points.len();
    let adj=build_adj(mesh,n);
    let data:Vec<f64>=adj.iter().map(|a|a.len() as f64).collect();
    let mut result=mesh.clone();
    result.point_data_mut().add_array(AnyDataArray::F64(DataArray::from_vec("Valence",data,1)));
    result
}

/// Extract n-ring as a separate PolyData (faces within the ring).
pub fn extract_n_ring_mesh(mesh: &PolyData, vertex: usize, radius: usize) -> PolyData {
    let ring_verts:std::collections::HashSet<usize>=n_ring(mesh,vertex,radius).into_iter().collect();
    let all_cells:Vec<Vec<i64>>=mesh.polys.iter().map(|c|c.to_vec()).collect();

    let mut pts=Points::<f64>::new();
    let mut polys=vtk_data::CellArray::new();
    let mut pm:std::collections::HashMap<usize,usize>=std::collections::HashMap::new();

    for cell in &all_cells{
        if cell.iter().all(|&pid|ring_verts.contains(&(pid as usize))){
            let mut ids=Vec::new();
            for &pid in cell{let old=pid as usize;
                let idx=*pm.entry(old).or_insert_with(||{let i=pts.len();pts.push(mesh.points.get(old));i});
                ids.push(idx as i64);
            }
            polys.push_cell(&ids);
        }
    }

    let mut result=PolyData::new();result.points=pts;result.polys=polys;result
}

/// Identify irregular vertices (valence != 6 for triangle meshes).
pub fn find_irregular_vertices(mesh: &PolyData) -> Vec<usize> {
    let n=mesh.points.len();
    let adj=build_adj(mesh,n);
    let boundary=find_boundary(mesh,n);
    (0..n).filter(|&i|{
        if adj[i].is_empty(){return false;}
        if boundary[i]{return adj[i].len()!=4;} // boundary regular = 4
        adj[i].len()!=6 // interior regular = 6
    }).collect()
}

fn build_adj(m:&PolyData,n:usize)->Vec<Vec<usize>>{
    let mut a:Vec<std::collections::HashSet<usize>>=vec![std::collections::HashSet::new();n];
    for c in m.polys.iter(){let nc=c.len();for i in 0..nc{
        let x=c[i] as usize;let y=c[(i+1)%nc] as usize;if x<n&&y<n{a[x].insert(y);a[y].insert(x);}
    }}a.into_iter().map(|s|s.into_iter().collect()).collect()
}

fn find_boundary(m:&PolyData,n:usize)->Vec<bool>{
    let mut ec:std::collections::HashMap<(usize,usize),usize>=std::collections::HashMap::new();
    for c in m.polys.iter(){let nc=c.len();for i in 0..nc{
        let a=c[i] as usize;let b=c[(i+1)%nc] as usize;*ec.entry((a.min(b),a.max(b))).or_insert(0)+=1;
    }}
    let mut bnd=vec![false;n];
    for (&(a,b),&c) in &ec{if c==1{bnd[a]=true;bnd[b]=true;}}
    bnd
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn one_ring_test() {
        let mesh=PolyData::from_triangles(
            vec![[0.0,0.0,0.0],[1.0,0.0,0.0],[0.5,1.0,0.0],[1.5,1.0,0.0]],
            vec![[0,1,2],[1,3,2]]);
        let ring=one_ring(&mesh,1);
        assert!(ring.len()>=2); // connected to at least 0,2,3
    }
    #[test]
    fn n_ring_test() {
        let mut pts=Vec::new();let mut tris=Vec::new();
        for y in 0..5{for x in 0..5{pts.push([x as f64,y as f64,0.0]);}}
        for y in 0..4{for x in 0..4{let bl=y*5+x;tris.push([bl,bl+1,bl+6]);tris.push([bl,bl+6,bl+5]);}}
        let mesh=PolyData::from_triangles(pts,tris);
        let ring=n_ring(&mesh,12,2); // center vertex, 2-ring
        assert!(ring.len()>4);
    }
    #[test]
    fn valence() {
        let mesh=PolyData::from_triangles(
            vec![[0.0,0.0,0.0],[1.0,0.0,0.0],[0.5,1.0,0.0]],vec![[0,1,2]]);
        let result=vertex_valence(&mesh);
        let arr=result.point_data().get_array("Valence").unwrap();
        let mut buf=[0.0f64]; arr.tuple_as_f64(0,&mut buf);
        assert_eq!(buf[0],2.0); // each vertex connects to 2 others
    }
    #[test]
    fn extract_ring() {
        let mut pts=Vec::new();let mut tris=Vec::new();
        for y in 0..5{for x in 0..5{pts.push([x as f64,y as f64,0.0]);}}
        for y in 0..4{for x in 0..4{let bl=y*5+x;tris.push([bl,bl+1,bl+6]);tris.push([bl,bl+6,bl+5]);}}
        let mesh=PolyData::from_triangles(pts,tris);
        let ring_mesh=extract_n_ring_mesh(&mesh,12,1);
        assert!(ring_mesh.polys.num_cells()>0);
        assert!(ring_mesh.polys.num_cells()<mesh.polys.num_cells());
    }
}
