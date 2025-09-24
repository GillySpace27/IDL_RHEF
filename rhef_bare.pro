pro rhef_bare, f
  compile_opt idl2
  read_sdo, f, hdr, im
  t = systime(1)
  x = findgen(hdr.naxis1) - hdr.crpix1
  y = findgen(hdr.naxis2) - hdr.crpix2
  xx = rebin(x, hdr.naxis1, hdr.naxis2)
  yy = rebin(transpose(y), hdr.naxis1, hdr.naxis2)
  r = sqrt(xx ^ 2 + yy ^ 2)
  h = histogram(r, binsize = 1, rev = revind, loc = rh)
  nh = n_elements(h)
  imout = im * 0.
  for index_hist = 0, nh - 1 do begin
    if h[index_hist] eq 0 then continue
    res = revind[revind[index_hist] : revind[index_hist + 1] - 1]
    n = n_elements(res)
    imnow = im[res]
    imout[res] = interpol(findgen(n) / (n - 1), imnow[sort(imnow)], imnow)
  endfor
  yl = 0.7
  yh = 0.4
  indl = where(imout le 0.5, comp = indh)
  imout[indl] = 0.5 * ((2 * imout[indl]) ^ yl)
  imout[indh] = 0.5 * (2 - (2 * (1 - imout[indh])) ^ yh)
  print, 'Time elapsed: ', systime(1) - t
  tvscl, rebin(imout, 1024, 1024)
end
