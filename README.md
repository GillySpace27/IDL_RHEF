# IDL_RHEF
An IDL implementation of the Radial Histogram Equalizing Filter (RHEF)

```idl
;------------------------------------------------------------------------------
; @file        rhef.pro
; @author      Chris Gilly + Anonymous Reviewer
; @brief       Radial Histogram Equalizing Filter (RHEF).
;              Reads an image, computes radius per pixel, performs a
;              per-annulus histogram equalization (monotonic remap) to
;              enhance coronal structure, then applies a two-part tonal
;              curve and displays a quicklook.
;
; @usage       IDL> .run rhef
;              IDL> rhef, f
;
; @param f     [in]  Filename or file-specifier accepted by READ_SDO.
;                    Must point to a single-image FITS (or similar) file.
;
; @keywords YL [in]  Optional. Low-half gamma exponent (float). Default = 0.7
;           YH [in]  Optional. High-half gamma exponent (float). Default = 0.4
;
; @notes       - This is a prototype for scientific distribution (GitHub).
;              - Requires READ_SDO to populate HDR (header) and IM (image).
;              - Uses IDL2 semantics for safer array handling (COMPILE_OPT IDL2).
;              - Quicklook display uses TVSCL with simple 1024x1024 resize.
;
; @example     IDL> rhef, 'aia.lev1_euv_12s_2024-01-01T00_00_00Z_193.fits', YL=0.6, YH=0.45
;
; @copyright   2025. Distribution under the MIT license.
;              Please cite appropriately if used in publications.
;------------------------------------------------------------------------------
```