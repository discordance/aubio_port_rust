use ::libc;
extern "C" {
    pub type _aubio_filterbank_t;
    pub type _aubio_dct_t;
    #[no_mangle]
    fn calloc(_: libc::c_ulong, _: libc::c_ulong) -> *mut libc::c_void;
    #[no_mangle]
    fn free(_: *mut libc::c_void);
    /* file interface */
    /* strings */
    /* Error reporting */
    /* Logging */
    /* * internal logging function, defined in utils/log.c */
    #[no_mangle]
    fn aubio_log(level: sint_t, fmt: *const char_t, _: ...) -> uint_t;
    #[no_mangle]
    fn fvec_mul(v: *mut fvec_t, s: smpl_t);
    /* * fvec_t buffer creation function

  \param length the length of the buffer to create

*/
    #[no_mangle]
    fn new_fvec(length: uint_t) -> *mut fvec_t;
    /* * fvec_t buffer deletion function

  \param s buffer to delete as returned by new_fvec()

*/
    #[no_mangle]
    fn del_fvec(s: *mut fvec_t);
    /* * make a copy of a vector

  \param s source vector
  \param t vector to copy to

*/
    #[no_mangle]
    fn fvec_copy(s: *const fvec_t, t: *mut fvec_t);
    /* * compute the \f$log10(x)\f$ of each vector elements

  \param s vector to modify

*/
    #[no_mangle]
    fn fvec_log10(s: *mut fvec_t);
    /* * create filterbank object

  \param n_filters number of filters to create
  \param win_s size of analysis buffer (and length the FFT transform)

*/
    #[no_mangle]
    fn new_aubio_filterbank(n_filters: uint_t, win_s: uint_t)
     -> *mut aubio_filterbank_t;
    /* * destroy filterbank object

  \param f filterbank object, as returned by new_aubio_filterbank()

*/
    #[no_mangle]
    fn del_aubio_filterbank(f: *mut aubio_filterbank_t);
    /* * compute filterbank

  \param f filterbank object, as returned by new_aubio_filterbank()
  \param in input spectrum containing an input spectrum of length `win_s`
  \param out output vector containing the energy found in each band, `nfilt` output values

*/
    #[no_mangle]
    fn aubio_filterbank_do(f: *mut aubio_filterbank_t, in_0: *const cvec_t,
                           out: *mut fvec_t);
    /* * set power parameter

  \param f filterbank object, as returned by new_aubio_filterbank()
  \param power Raise norm of the input spectrum norm to this power before
  computing filterbank.  Defaults to `1`.

 */
    #[no_mangle]
    fn aubio_filterbank_set_power(f: *mut aubio_filterbank_t, power: smpl_t)
     -> uint_t;
    /* * get power parameter

  \param f filterbank object, as returned by new_aubio_filterbank()
  \return current power parameter. Defaults to `1`.

 */
    #[no_mangle]
    fn aubio_filterbank_get_power(f: *mut aubio_filterbank_t) -> smpl_t;
    /* * filterbank initialization for Mel filters using Slaney's coefficients

  \param fb filterbank object
  \param samplerate audio sampling rate, in Hz

  The filter coefficients are built to match exactly Malcolm Slaney's Auditory
  Toolbox implementation (see file mfcc.m). The number of filters should be 40.

  References
  ----------

  Malcolm Slaney, *Auditory Toolbox Version 2, Technical Report #1998-010*
  https://engineering.purdue.edu/~malcolm/interval/1998-010/

*/
    #[no_mangle]
    fn aubio_filterbank_set_mel_coeffs_slaney(fb: *mut aubio_filterbank_t,
                                              samplerate: smpl_t) -> uint_t;
    /* * Mel filterbank initialization

  \param fb filterbank object
  \param samplerate audio sampling rate
  \param fmin start frequency, in Hz
  \param fmax end frequency, in Hz

  The filterbank will be initialized with bands linearly spaced in the mel
  scale, from `fmin` to `fmax`.

  References
  ----------

  Malcolm Slaney, *Auditory Toolbox Version 2, Technical Report #1998-010*
  https://engineering.purdue.edu/~malcolm/interval/1998-010/

*/
    #[no_mangle]
    fn aubio_filterbank_set_mel_coeffs(fb: *mut aubio_filterbank_t,
                                       samplerate: smpl_t, fmin: smpl_t,
                                       fmax: smpl_t) -> uint_t;
    /* * Mel filterbank initialization

  \param fb filterbank object
  \param samplerate audio sampling rate
  \param fmin start frequency, in Hz
  \param fmax end frequency, in Hz

  The bank of filters will be initalized to to cover linearly spaced bands in
  the Htk mel scale, from `fmin` to `fmax`.

  References
  ----------

  Douglas O'Shaughnessy (1987). *Speech communication: human and machine*.
  Addison-Wesley. p. 150. ISBN 978-0-201-16520-3.

  HTK Speech Recognition Toolkit: http://htk.eng.cam.ac.uk/

*/
    #[no_mangle]
    fn aubio_filterbank_set_mel_coeffs_htk(fb: *mut aubio_filterbank_t,
                                           samplerate: smpl_t, fmin: smpl_t,
                                           fmax: smpl_t) -> uint_t;
    /* * create new DCT computation object

  \param size length of the DCT

*/
    #[no_mangle]
    fn new_aubio_dct(size: uint_t) -> *mut aubio_dct_t;
    /* * compute forward DCT

  \param s dct object as returned by new_aubio_dct
  \param input input signal
  \param dct_output transformed input array

*/
    #[no_mangle]
    fn aubio_dct_do(s: *mut aubio_dct_t, input: *const fvec_t,
                    dct_output: *mut fvec_t);
    /* * delete DCT object

  \param s dct object as returned by new_aubio_dct

*/
    #[no_mangle]
    fn del_aubio_dct(s: *mut aubio_dct_t);
}
/*
  Copyright (C) 2003-2015 Paul Brossier <piem@aubio.org>

  This file is part of aubio.

  aubio is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  aubio is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with aubio.  If not, see <http://www.gnu.org/licenses/>.

*/
/* * \file

  Definition of data types used in aubio

*/
/* * defined to 1 if aubio is compiled in double precision */
/* * short sample format (32 or 64 bits) */
pub type smpl_t = libc::c_float;
/* * print format for sample in double precision */
/* * unsigned integer */
pub type uint_t = libc::c_uint;
/* * signed integer */
pub type sint_t = libc::c_int;
/* * character */
pub type char_t = libc::c_char;
/*
  Copyright (C) 2003-2015 Paul Brossier <piem@aubio.org>

  This file is part of aubio.

  aubio is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  aubio is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with aubio.  If not, see <http://www.gnu.org/licenses/>.

*/
/* * \file

  Vector of real-valued data

  This file specifies the ::fvec_t buffer type, which is used throughout aubio
  to store vector of real-valued ::smpl_t.

  \example test-fvec.c

*/
/* * Buffer for real data

  Vector of real-valued data

  ::fvec_t is is the structure used to store vector of real-valued data, ::smpl_t .

  \code

  uint_t buffer_size = 1024;

  // create a vector of 512 values
  fvec_t * input = new_fvec (buffer_size);

  // set some values of the vector
  input->data[23] = 2.;
  // ..

  // compute the mean of the vector
  mean = fvec_mean(a_vector);

  // destroy the vector
  del_fvec(a_vector);

  \endcode

  See `examples/` and `tests/src` directories for more examples.

 */
#[derive(Copy, Clone)]
#[repr(C)]
pub struct fvec_t {
    pub length: uint_t,
    pub data: *mut smpl_t,
}
pub type C2RustUnnamed = libc::c_uint;
pub const AUBIO_FAIL: C2RustUnnamed = 1;
pub const AUBIO_OK: C2RustUnnamed = 0;
pub type aubio_log_level = libc::c_uint;
pub const AUBIO_LOG_LAST_LEVEL: aubio_log_level = 5;
pub const AUBIO_LOG_WRN: aubio_log_level = 4;
pub const AUBIO_LOG_DBG: aubio_log_level = 3;
pub const AUBIO_LOG_MSG: aubio_log_level = 2;
pub const AUBIO_LOG_INF: aubio_log_level = 1;
pub const AUBIO_LOG_ERR: aubio_log_level = 0;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct cvec_t {
    pub length: uint_t,
    pub norm: *mut smpl_t,
    pub phas: *mut smpl_t,
}
/*
  Copyright (C) 2007-2013 Paul Brossier <piem@aubio.org>
                      and Amaury Hazan <ahazan@iua.upf.edu>

  This file is part of aubio.

  aubio is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  aubio is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with aubio.  If not, see <http://www.gnu.org/licenses/>.

*/
/* * \file

  Filterbank object

  General-purpose spectral filterbank object.

  \example spectral/test-filterbank.c

*/
/* * filterbank object

  This object stores a matrix of spectral filter coefficients.

 */
pub type aubio_filterbank_t = _aubio_filterbank_t;
/*
  Copyright (C) 2017 Paul Brossier <piem@aubio.org>

  This file is part of aubio.

  aubio is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  aubio is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with aubio.  If not, see <http://www.gnu.org/licenses/>.

*/
/* * \file

  Discrete Cosine Transform

  Functions aubio_dct_do() and aubio_dct_rdo() are equivalent to MATLAB/Octave
  dct() and idct() functions, as well as scipy.fftpack.dct(x, norm='ortho') and
  scipy.fftpack.idct(x, norm='ortho')

  \example spectral/test-dct.c

*/
/* * DCT object

  This object computes forward and backward DCT type 2 with orthonormal
  scaling.

*/
pub type aubio_dct_t = _aubio_dct_t;
/*
  Copyright (C) 2007-2009 Paul Brossier <piem@aubio.org>
                      and Amaury Hazan <ahazan@iua.upf.edu>

  This file is part of aubio.

  aubio is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  aubio is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with aubio.  If not, see <http://www.gnu.org/licenses/>.

*/
/* * Internal structure for mfcc object */
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _aubio_mfcc_t {
    pub win_s: uint_t,
    pub samplerate: uint_t,
    pub n_filters: uint_t,
    pub n_coefs: uint_t,
    pub fb: *mut aubio_filterbank_t,
    pub in_dct: *mut fvec_t,
    pub dct: *mut aubio_dct_t,
    pub output: *mut fvec_t,
    pub scale: smpl_t,
}
/*
  Copyright (C) 2007-2013 Paul Brossier <piem@aubio.org>
                      and Amaury Hazan <ahazan@iua.upf.edu>

  This file is part of aubio.

  aubio is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  aubio is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with aubio.  If not, see <http://www.gnu.org/licenses/>.

*/
/* * \file

  Mel-Frequency Cepstrum Coefficients object

  This object computes MFCC coefficients on an input cvec_t.

  The implementation follows the specifications established by Malcolm Slaney
  in its Auditory Toolbox, available online at the following address (see
  file mfcc.m):

  https://engineering.purdue.edu/~malcolm/interval/1998-010/

  \example spectral/test-mfcc.c

*/
/* * mfcc object */
pub type aubio_mfcc_t = _aubio_mfcc_t;
/* * create mfcc object

  \param buf_size size of analysis buffer (and length the FFT transform)
  \param samplerate audio sampling rate
  \param n_coeffs number of desired coefficients
  \param n_filters number of desired filters

*/
#[no_mangle]
pub unsafe extern "C" fn new_aubio_mfcc(mut win_s: uint_t,
                                        mut n_filters: uint_t,
                                        mut n_coefs: uint_t,
                                        mut samplerate: uint_t)
 -> *mut aubio_mfcc_t {
    /* allocate space for mfcc object */
    let mut mfcc: *mut aubio_mfcc_t =
        calloc(::std::mem::size_of::<aubio_mfcc_t>() as libc::c_ulong,
               1 as libc::c_int as libc::c_ulong) as *mut aubio_mfcc_t;
    if n_coefs as sint_t <= 0 as libc::c_int {
        aubio_log(AUBIO_LOG_ERR as libc::c_int,
                  b"AUBIO ERROR: mfcc: n_coefs should be > 0, got %d\n\x00" as
                      *const u8 as *const libc::c_char, n_coefs);
    } else if samplerate as sint_t <= 0 as libc::c_int {
        aubio_log(AUBIO_LOG_ERR as libc::c_int,
                  b"AUBIO ERROR: mfcc: samplerate should be > 0, got %d\n\x00"
                      as *const u8 as *const libc::c_char, samplerate);
    } else {
        (*mfcc).win_s = win_s;
        (*mfcc).samplerate = samplerate;
        (*mfcc).n_filters = n_filters;
        (*mfcc).n_coefs = n_coefs;
        /* filterbank allocation */
        (*mfcc).fb = new_aubio_filterbank(n_filters, (*mfcc).win_s);
        if !(*mfcc).fb.is_null() {
            if n_filters == 40 as libc::c_int as libc::c_uint {
                aubio_filterbank_set_mel_coeffs_slaney((*mfcc).fb,
                                                       samplerate as smpl_t);
            } else {
                aubio_filterbank_set_mel_coeffs((*mfcc).fb,
                                                samplerate as smpl_t,
                                                0 as libc::c_int as smpl_t,
                                                (samplerate as libc::c_double
                                                     / 2.0f64) as smpl_t);
            }
            /* allocating buffers */
            (*mfcc).in_dct = new_fvec(n_filters);
            (*mfcc).dct = new_aubio_dct(n_filters);
            (*mfcc).output = new_fvec(n_filters);
            if !((*mfcc).in_dct.is_null() || (*mfcc).dct.is_null() ||
                     (*mfcc).output.is_null()) {
                (*mfcc).scale = 1.0f64 as smpl_t;
                return mfcc
            }
        }
    }
    del_aubio_mfcc(mfcc);
    return 0 as *mut aubio_mfcc_t;
}
/* * delete mfcc object

  \param mf mfcc object as returned by new_aubio_mfcc

*/
#[no_mangle]
pub unsafe extern "C" fn del_aubio_mfcc(mut mf: *mut aubio_mfcc_t) {
    if !(*mf).fb.is_null() { del_aubio_filterbank((*mf).fb); }
    if !(*mf).in_dct.is_null() { del_fvec((*mf).in_dct); }
    if !(*mf).dct.is_null() { del_aubio_dct((*mf).dct); }
    if !(*mf).output.is_null() { del_fvec((*mf).output); }
    free(mf as *mut libc::c_void);
}
/* * mfcc object processing

  \param mf mfcc object as returned by new_aubio_mfcc
  \param in input spectrum (buf_size long)
  \param out output mel coefficients buffer (n_coeffs long)

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_mfcc_do(mut mf: *mut aubio_mfcc_t,
                                       mut in_0: *const cvec_t,
                                       mut out: *mut fvec_t) {
    let mut tmp: fvec_t = fvec_t{length: 0, data: 0 as *mut smpl_t,};
    /* compute filterbank */
    aubio_filterbank_do((*mf).fb, in_0, (*mf).in_dct);
    /* compute log10 */
    fvec_log10((*mf).in_dct);
    if (*mf).scale != 1 as libc::c_int as libc::c_float {
        fvec_mul((*mf).in_dct, (*mf).scale);
    }
    /* compute mfccs */
    aubio_dct_do((*mf).dct, (*mf).in_dct, (*mf).output);
    // copy only first n_coeffs elements
  // TODO assert mf->output->length == n_coeffs
    tmp.data = (*(*mf).output).data;
    tmp.length = (*out).length;
    fvec_copy(&mut tmp, out);
}
/* * set power parameter

  \param mf mfcc object, as returned by new_aubio_mfcc()
  \param power Raise norm of the input spectrum norm to this power before
  computing filterbank.  Defaults to `1`.

  See aubio_filterbank_set_power().

 */
#[no_mangle]
pub unsafe extern "C" fn aubio_mfcc_set_power(mut mf: *mut aubio_mfcc_t,
                                              mut power: smpl_t) -> uint_t {
    return aubio_filterbank_set_power((*mf).fb, power);
}
/* * get power parameter

  \param mf mfcc object, as returned by new_aubio_mfcc()
  \return current power parameter. Defaults to `1`.

  See aubio_filterbank_get_power().

 */
#[no_mangle]
pub unsafe extern "C" fn aubio_mfcc_get_power(mut mf: *mut aubio_mfcc_t)
 -> smpl_t {
    return aubio_filterbank_get_power((*mf).fb);
}
/* * set scaling parameter

  \param mf mfcc object, as returned by new_aubio_mfcc()
  \param scale Scaling value to apply.

  Scales the output of the filterbank after taking its logarithm and before
  computing the DCT. Defaults to `1`.

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_mfcc_set_scale(mut mf: *mut aubio_mfcc_t,
                                              mut scale: smpl_t) -> uint_t {
    (*mf).scale = scale;
    return AUBIO_OK as libc::c_int as uint_t;
}
/* * get scaling parameter

  \param mf mfcc object, as returned by new_aubio_mfcc()
  \return current scaling parameter. Defaults to `1`.

 */
#[no_mangle]
pub unsafe extern "C" fn aubio_mfcc_get_scale(mut mf: *mut aubio_mfcc_t)
 -> smpl_t {
    return (*mf).scale;
}
/* * Mel filterbank initialization

  \param mf mfcc object
  \param fmin start frequency, in Hz
  \param fmax end frequency, in Hz

  The filterbank will be initialized with bands linearly spaced in the mel
  scale, from `fmin` to `fmax`.

  See also
  --------

  aubio_filterbank_set_mel_coeffs()

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_mfcc_set_mel_coeffs(mut mf: *mut aubio_mfcc_t,
                                                   mut freq_min: smpl_t,
                                                   mut freq_max: smpl_t)
 -> uint_t {
    return aubio_filterbank_set_mel_coeffs((*mf).fb,
                                           (*mf).samplerate as smpl_t,
                                           freq_min, freq_max);
}
/* * Mel filterbank initialization

  \param mf mfcc object
  \param fmin start frequency, in Hz
  \param fmax end frequency, in Hz

  The bank of filters will be initalized to to cover linearly spaced bands in
  the Htk mel scale, from `fmin` to `fmax`.

  See also
  --------

  aubio_filterbank_set_mel_coeffs_htk()

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_mfcc_set_mel_coeffs_htk(mut mf:
                                                           *mut aubio_mfcc_t,
                                                       mut freq_min: smpl_t,
                                                       mut freq_max: smpl_t)
 -> uint_t {
    return aubio_filterbank_set_mel_coeffs_htk((*mf).fb,
                                               (*mf).samplerate as smpl_t,
                                               freq_min, freq_max);
}
/* * Mel filterbank initialization (Auditory Toolbox's parameters)

  \param mf mfcc object

  The filter coefficients are built to match exactly Malcolm Slaney's Auditory
  Toolbox implementation. The number of filters should be 40.

  This is the default filterbank when `mf` was created with `n_filters = 40`.

  See also
  --------

  aubio_filterbank_set_mel_coeffs_slaney()

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_mfcc_set_mel_coeffs_slaney(mut mf:
                                                              *mut aubio_mfcc_t)
 -> uint_t {
    return aubio_filterbank_set_mel_coeffs_slaney((*mf).fb,
                                                  (*mf).samplerate as smpl_t);
}
