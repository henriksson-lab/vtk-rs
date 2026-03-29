//! Compute percentile values of scalar arrays.
use vtk_data::{AnyDataArray, DataArray, PolyData};
pub fn scalar_percentile(mesh: &PolyData, array_name: &str, percentile: f64) -> f64 {
    let arr=match mesh.point_data().get_array(array_name){Some(a) if a.num_components()==1=>a,_=>return 0.0};
    let n=arr.num_tuples();if n==0{return 0.0;}let mut buf=[0.0f64];
    let mut vals:Vec<f64>=(0..n).map(|i|{arr.tuple_as_f64(i,&mut buf);buf[0]}).collect();
    vals.sort_by(|a,b|a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let idx=((percentile/100.0*n as f64).floor() as usize).min(n-1);vals[idx]
}
pub fn clip_by_percentile(mesh: &PolyData, array_name: &str, lo_pct: f64, hi_pct: f64) -> PolyData {
    let lo=scalar_percentile(mesh,array_name,lo_pct);
    let hi=scalar_percentile(mesh,array_name,hi_pct);
    let arr=match mesh.point_data().get_array(array_name){Some(a) if a.num_components()==1=>a,_=>return mesh.clone()};
    let n=arr.num_tuples();let mut buf=[0.0f64];
    let data:Vec<f64>=(0..n).map(|i|{arr.tuple_as_f64(i,&mut buf);buf[0].clamp(lo,hi)}).collect();
    let mut r=mesh.clone();
    r.point_data_mut().add_array(AnyDataArray::F64(DataArray::from_vec(array_name,data,1)));r
}
pub fn normalize_by_percentile(mesh: &PolyData, array_name: &str, lo_pct: f64, hi_pct: f64) -> PolyData {
    let lo=scalar_percentile(mesh,array_name,lo_pct);
    let hi=scalar_percentile(mesh,array_name,hi_pct);
    let range=(hi-lo).max(1e-15);
    let arr=match mesh.point_data().get_array(array_name){Some(a) if a.num_components()==1=>a,_=>return mesh.clone()};
    let n=arr.num_tuples();let mut buf=[0.0f64];
    let data:Vec<f64>=(0..n).map(|i|{arr.tuple_as_f64(i,&mut buf);((buf[0]-lo)/range).clamp(0.0,1.0)}).collect();
    let mut r=mesh.clone();
    r.point_data_mut().add_array(AnyDataArray::F64(DataArray::from_vec(array_name,data,1)));r
}
#[cfg(test)] mod tests { use super::*;
    #[test] fn test_pct() { let mut m=PolyData::from_triangles(vec![[0.0,0.0,0.0],[1.0,0.0,0.0],[0.5,1.0,0.0]],vec![[0,1,2]]);
        m.point_data_mut().add_array(AnyDataArray::F64(DataArray::from_vec("s",vec![0.0,50.0,100.0],1)));
        let p50=scalar_percentile(&m,"s",50.0); assert!((p50-50.0).abs()<1e-10); }
    #[test] fn test_clip() { let mut m=PolyData::from_triangles(vec![[0.0,0.0,0.0],[1.0,0.0,0.0],[0.5,1.0,0.0]],vec![[0,1,2]]);
        m.point_data_mut().add_array(AnyDataArray::F64(DataArray::from_vec("s",vec![0.0,50.0,100.0],1)));
        let r=clip_by_percentile(&m,"s",10.0,90.0); assert!(r.point_data().get_array("s").is_some()); }
    #[test] fn test_norm() { let mut m=PolyData::from_triangles(vec![[0.0,0.0,0.0],[1.0,0.0,0.0],[0.5,1.0,0.0]],vec![[0,1,2]]);
        m.point_data_mut().add_array(AnyDataArray::F64(DataArray::from_vec("s",vec![10.0,50.0,90.0],1)));
        let r=normalize_by_percentile(&m,"s",0.0,100.0); assert!(r.point_data().get_array("s").is_some()); } }
