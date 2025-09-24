;------------------------------------------------------------------------------
; @file        rhef.pro
; @author      Gilly Gilbert
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
; @copyright   2025. Distribution under an open-science license of your
;              choosing. Please cite appropriately if used in publications.
;------------------------------------------------------------------------------
pro rhef, f, YL=yl, YH=yh
    compile_opt idl2

    ;-----------------------------
    ; Validate inputs & set defaults
    ;-----------------------------
    if n_params() lt 1 then begin
        message, 'Usage: RHEF, f [, YL=yl, YH=yh]  ; f = filename accepted by READ_SDO', /info
        return
    endif

    ; Default exponents if not supplied by caller
    if n_elements(yl) eq 0 then yl = 0.7
    if n_elements(yh) eq 0 then yh = 0.4

    ;-----------------------------
    ; Ingest data (header + image)
    ;-----------------------------
    read_sdo, f, hdr, im

    ; Simple sanity checks
    if ~size(/type, im) then message, 'RHEF: IM is undefined after READ_SDO.'
    if ~keyword_set(hdr) then message, 'RHEF: HDR is undefined after READ_SDO.'

    ;-----------------------------
    ; Start a timer for basic performance feedback
    ;-----------------------------
    t = systime(1)

    ;-----------------------------
    ; Build pixel-centered coordinate grids (x, y) in pixels relative to CRPIX
    ;   x runs along naxis1, y along naxis2; then broadcast to full image size.
    ;-----------------------------
    x = findgen(hdr.naxis1) - hdr.crpix1
    y = findgen(hdr.naxis2) - hdr.crpix2
    xx = rebin(x,              hdr.naxis1, hdr.naxis2)           ; broadcast x
    yy = rebin(transpose(y),   hdr.naxis1, hdr.naxis2)           ; broadcast y

    ;-----------------------------
    ; Radial distance per pixel (in pixel units)
    ;-----------------------------
    r = sqrt(xx^2 + yy^2)

    ;-----------------------------
    ; Build a 1-pixel-wide radial histogram with reverse indices
    ;   - RH: bin centers (locations)
    ;   - H:  counts per bin
    ;   - REVIND: reverse indices to map bins -> pixel indices
    ;-----------------------------
    h = histogram(r, binsize=1, rev=revind, loc=rh)
    nh = n_elements(h)

    ;-----------------------------
    ; Allocate output image (same type as IM) and initialize to zero
    ;-----------------------------
    imout = im * 0.

    ;-----------------------------
    ; Per-annulus monotonic remap (histogram equalization within each ring)
    ;   For each radial bin:
    ;     * Gather pixel indices (RES)
    ;     * Sort intensities and map to [0,1] via rank / (N-1)
    ;     * Interpolate back to original ordering so output preserves topology
    ;-----------------------------
    for index_hist = 0, nh - 1 do begin
        if h[index_hist] eq 0 then continue  ; empty ring → skip

        ; Reverse indices range for this ring
        res = revind[ revind[index_hist] : (revind[index_hist+1] - 1) ]
        n   = n_elements(res)

        imnow = im[res]
        ; Create a monotonically increasing mapping from sorted values to [0,1]
        ; NOTE: findgen(n)/(n-1) maps ranks to [0,1]; INTERPOL applies it back.
        imout[res] = interpol( findgen(n)/(n-1), imnow[ sort(imnow) ], imnow )
    endfor

    ;-----------------------------
    ; Two-part tonal curve: compress shadows & highlights with exponents
    ;   controlled by YL (low) and YH (high). Work on [0,1] split at 0.5.
    ;-----------------------------
    indl = where(imout le 0.5, comp = indh)   ; INDL: low half, INDH: complement

    ; Low half remap: scale to [0,1], apply exponent YL, then rescale back
    imout[indl] = 0.5 * ((2 * imout[indl])^yl)

    ; High half remap: symmetric treatment using exponent YH
    imout[indh] = 0.5 * ( 2 - ( 2 * (1 - imout[indh]) )^yh )

    ;-----------------------------
    ; Report runtime
    ;-----------------------------
    print, 'RHEF: Time elapsed (s): ', systime(1) - t

    ;-----------------------------
    ; Quicklook display (scaled to 1024x1024 for convenience)
    ;-----------------------------
    tvscl, rebin(imout, 1024, 1024)

end