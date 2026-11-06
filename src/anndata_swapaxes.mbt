/// anndata_swapaxes.mbt — Swap axes of AnnData
/// Item 11: swapaxes

///|
/// Swap the axes of an AnnData. For a 2D AnnData, this is equivalent to
/// transpose but more general (supports N-D in principle).
/// axis1 and axis2 are 0-based axis indices (0=obs, 1=var).
pub fn anndata_swapaxes(
  adata : AnnData,
  axis1~ : Int32 = 0,
  axis2~ : Int32 = 1,
) -> AnnData {
  // For 2D AnnData, swapping axes 0 and 1 is transpose
  if (axis1 == 0 && axis2 == 1) || (axis1 == 1 && axis2 == 0) {
    anndata_transpose(adata)
  } else {
    // Same axis or invalid → return copy
    adata.copy()
  }
}

///|
/// Swap obs and var annotations (more thorough than transpose:
/// also swaps obsm↔varm, obsp↔varp).
pub fn anndata_swap_axes_full(adata : AnnData) -> AnnData {
  let result = anndata_transpose(adata)
  // Swap obsm ↔ varm
  let obsm_keys = adata.obsm.keys().collect()
  for k in obsm_keys {
    result.varm[k] = adata.obsm[k]
  }
  let varm_keys = adata.varm.keys().collect()
  for k in varm_keys {
    result.obsm[k] = adata.varm[k]
  }
  // Swap obsp ↔ varp
  let obsp_keys = adata.obsp.keys().collect()
  for k in obsp_keys {
    result.varp[k] = adata.obsp[k]
  }
  let varp_keys = adata.varp.keys().collect()
  for k in varp_keys {
    result.obsp[k] = adata.varp[k]
  }
  result
}
