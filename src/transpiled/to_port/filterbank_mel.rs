use ::libc;
extern "C" {
    pub type _aubio_filterbank_t;
    #[no_mangle]
    fn powf(_: libc::c_float, _: libc::c_float) -> libc::c_float;
    #[no_mangle]
    fn aubio_log(level: sint_t, fmt: *const char_t, _: ...) -> uint_t;
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
    /* * set all elements to ones

  \param s vector to modify

*/
    #[no_mangle]
    fn fvec_ones(s: *mut fvec_t);
    /* * convert frequency bin to frequency (Hz) */
    #[no_mangle]
    fn aubio_bintofreq(bin: smpl_t, samplerate: smpl_t, fftsize: smpl_t)
     -> smpl_t;
    /* * convert frequency (Hz) to mel

  \param freq input frequency, in Hz

  \return output mel

  Converts a scalar from the frequency domain to the mel scale using Slaney
  Auditory Toolbox's implementation:

  If \f$ f < 1000 \f$, \f$ m = 3 f / 200 \f$.

  If \f$ f >= 1000 \f$, \f$ m = 1000 + 27 \frac{{ln}(f) - ln(1000))}
  {{ln}(6400) - ln(1000)}
  \f$

  See also
  --------

  aubio_meltohz(), aubio_hztomel_htk().

*/
    #[no_mangle]
    fn aubio_hztomel(freq: smpl_t) -> smpl_t;
    /* * convert mel to frequency (Hz)

  \param mel input mel

  \return output frequency, in Hz

  Converts a scalar from the mel scale to the frequency domain using Slaney
  Auditory Toolbox's implementation:

  If \f$ f < 1000 \f$, \f$ f = 200 m/3 \f$.

  If \f$ f \geq 1000 \f$, \f$ f = 1000 + \left(\frac{6400}{1000}\right)
  ^{\frac{m - 1000}{27}} \f$

  See also
  --------

  aubio_hztomel(), aubio_meltohz_htk().

  References
  ----------

  Malcolm Slaney, *Auditory Toolbox Version 2, Technical Report #1998-010*
  https://engineering.purdue.edu/~malcolm/interval/1998-010/

*/
    #[no_mangle]
    fn aubio_meltohz(mel: smpl_t) -> smpl_t;
    /* * convert frequency (Hz) to mel

  \param freq input frequency, in Hz

  \return output mel

  Converts a scalar from the frequency domain to the mel scale, using the
  equation defined by O'Shaughnessy, as implemented in the HTK speech
  recognition toolkit:

  \f$ m = 1127 + ln(1 + \frac{f}{700}) \f$

  See also
  --------

  aubio_meltohz_htk(), aubio_hztomel().

  References
  ----------

  Douglas O'Shaughnessy (1987). *Speech communication: human and machine*.
  Addison-Wesley. p. 150. ISBN 978-0-201-16520-3.

  HTK Speech Recognition Toolkit: http://htk.eng.cam.ac.uk/

 */
    #[no_mangle]
    fn aubio_hztomel_htk(freq: smpl_t) -> smpl_t;
    /* * convert mel to frequency (Hz)

  \param mel input mel

  \return output frequency, in Hz

  Converts a scalar from the mel scale to the frequency domain, using the
  equation defined by O'Shaughnessy, as implemented in the HTK speech
  recognition toolkit:

  \f$ f = 700 * {e}^\left(\frac{f}{1127} - 1\right) \f$

  See also
  --------

  aubio_hztomel_htk(), aubio_meltohz().

*/
    #[no_mangle]
    fn aubio_meltohz_htk(mel: smpl_t) -> smpl_t;
    /* * set all elements to zero

  \param s vector to modify

*/
    #[no_mangle]
    fn fmat_zeros(s: *mut fmat_t);
    /* * return a pointer to the matrix object containing all filter coefficients

  \param f filterbank object, as returned by new_aubio_filterbank()

 */
    #[no_mangle]
    fn aubio_filterbank_get_coeffs(f: *const aubio_filterbank_t)
     -> *mut fmat_t;
    /* * get norm parameter

  \param f filterbank object, as returned by new_aubio_filterbank()
  \returns `1` if norm is set, `0` otherwise. Defaults to `1`.

 */
    #[no_mangle]
    fn aubio_filterbank_get_norm(f: *mut aubio_filterbank_t) -> smpl_t;
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
/*
  Copyright (C) 2016 Paul Brossier <piem@aubio.org>

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

  Logging features

  This file specifies ::aubio_log_set_function and
  ::aubio_log_set_level_function, which let you define one or several custom
  logging functions to redirect warnings and errors from aubio to your
  application. The custom function should have the prototype defined in
  ::aubio_log_function_t.

  After a call to ::aubio_log_set_level_function, ::aubio_log_reset can be used
  to reset each logging functions to the default ones.

  \example utils/test-log.c

*/
/* * list of logging levels */
pub type aubio_log_level = libc::c_uint;
/* *< number of valid levels */
/* *< warnings */
pub const AUBIO_LOG_LAST_LEVEL: aubio_log_level = 5;
/* *< debug messages */
pub const AUBIO_LOG_WRN: aubio_log_level = 4;
/* *< general messages */
pub const AUBIO_LOG_DBG: aubio_log_level = 3;
/* *< infos */
pub const AUBIO_LOG_MSG: aubio_log_level = 2;
/* *< critical errors */
pub const AUBIO_LOG_INF: aubio_log_level = 1;
pub const AUBIO_LOG_ERR: aubio_log_level = 0;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct fmat_t {
    pub length: uint_t,
    pub height: uint_t,
    pub data: *mut *mut smpl_t,
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

  Filterbank object coefficients initialization

  Functions to create set the ::aubio_filterbank_t coefficients to
    - ::aubio_filterbank_set_triangle_bands: overlapping triangular bands,
    - ::aubio_filterbank_set_mel_coeffs_slaney: Mel frequency bands.

  \example spectral/test-filterbank_mel.c

*/
/* * filterbank initialization with triangular and overlapping bands

  \param fb filterbank object
  \param freqs arbitrary array of boundary frequencies
  \param samplerate audio sampling rate

  This function computes the coefficients of the filterbank based on the
  boundaries found in freqs, in Hz, and using triangular overlapping bands.

*/
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
#[no_mangle]
pub unsafe extern "C" fn aubio_filterbank_set_triangle_bands(mut fb:
                                                                 *mut aubio_filterbank_t,
                                                             mut freqs:
                                                                 *const fvec_t,
                                                             mut samplerate:
                                                                 smpl_t)
 -> uint_t {
    let mut filters: *mut fmat_t =
        aubio_filterbank_get_coeffs(fb); /* filter counter */
    let mut n_filters: uint_t = (*filters).height; /* bin counter */
    let mut win_s: uint_t = (*filters).length;
    let mut lower_freqs: *mut fvec_t = 0 as *mut fvec_t;
    let mut upper_freqs: *mut fvec_t = 0 as *mut fvec_t;
    let mut center_freqs: *mut fvec_t = 0 as *mut fvec_t;
    let mut triangle_heights: *mut fvec_t = 0 as *mut fvec_t;
    let mut fft_freqs: *mut fvec_t = 0 as *mut fvec_t;
    let mut fn_0: uint_t = 0;
    let mut bin: uint_t = 0;
    let mut riseInc: smpl_t = 0.;
    let mut downInc: smpl_t = 0.;
    /* freqs define the bands of triangular overlapping windows.
     throw a warning if filterbank object fb is too short. */
    if (*freqs).length.wrapping_sub(2 as libc::c_int as libc::c_uint) >
           n_filters {
        aubio_log(AUBIO_LOG_WRN as libc::c_int,
                  b"AUBIO WARNING: not enough filters, %d allocated but %d requested\n\x00"
                      as *const u8 as *const libc::c_char, n_filters,
                  (*freqs).length.wrapping_sub(2 as libc::c_int as
                                                   libc::c_uint));
    }
    if (*freqs).length.wrapping_sub(2 as libc::c_int as libc::c_uint) <
           n_filters {
        aubio_log(AUBIO_LOG_WRN as libc::c_int,
                  b"AUBIO WARNING: too many filters, %d allocated but %d requested\n\x00"
                      as *const u8 as *const libc::c_char, n_filters,
                  (*freqs).length.wrapping_sub(2 as libc::c_int as
                                                   libc::c_uint));
    }
    fn_0 = 0 as libc::c_int as uint_t;
    while fn_0 < (*freqs).length {
        if *(*freqs).data.offset(fn_0 as isize) <
               0 as libc::c_int as libc::c_float {
            aubio_log(AUBIO_LOG_ERR as libc::c_int,
                      b"AUBIO ERROR: filterbank_mel: freqs must contain only positive values.\n\x00"
                          as *const u8 as *const libc::c_char);
            return AUBIO_FAIL as libc::c_int as uint_t
        } else {
            if *(*freqs).data.offset(fn_0 as isize) >
                   samplerate / 2 as libc::c_int as libc::c_float {
                aubio_log(AUBIO_LOG_WRN as libc::c_int,
                          b"AUBIO WARNING: filterbank_mel: freqs should contain only values < samplerate / 2.\n\x00"
                              as *const u8 as *const libc::c_char);
            } else if fn_0 > 0 as libc::c_int as libc::c_uint &&
                          *(*freqs).data.offset(fn_0 as isize) <
                              *(*freqs).data.offset(fn_0.wrapping_sub(1 as
                                                                          libc::c_int
                                                                          as
                                                                          libc::c_uint)
                                                        as isize) {
                aubio_log(AUBIO_LOG_ERR as libc::c_int,
                          b"AUBIO ERROR: filterbank_mel: freqs should be a list of frequencies sorted from low to high, but freq[%d] < freq[%d-1]\n\x00"
                              as *const u8 as *const libc::c_char, fn_0,
                          fn_0);
                return AUBIO_FAIL as libc::c_int as uint_t
            } else {
                if fn_0 > 0 as libc::c_int as libc::c_uint &&
                       *(*freqs).data.offset(fn_0 as isize) ==
                           *(*freqs).data.offset(fn_0.wrapping_sub(1 as
                                                                       libc::c_int
                                                                       as
                                                                       libc::c_uint)
                                                     as isize) {
                    aubio_log(AUBIO_LOG_WRN as libc::c_int,
                              b"AUBIO WARNING: filterbank_mel: set_triangle_bands received a list with twice the frequency %f\n\x00"
                                  as *const u8 as *const libc::c_char,
                              *(*freqs).data.offset(fn_0 as isize) as
                                  libc::c_double);
                }
            }
        }
        fn_0 = fn_0.wrapping_add(1)
    }
    /* convenience reference to lower/center/upper frequency for each triangle */
    lower_freqs = new_fvec(n_filters);
    upper_freqs = new_fvec(n_filters);
    center_freqs = new_fvec(n_filters);
    /* height of each triangle */
    triangle_heights = new_fvec(n_filters);
    /* lookup table of each bin frequency in hz */
    fft_freqs = new_fvec(win_s);
    /* fill up the lower/center/upper */
    fn_0 = 0 as libc::c_int as uint_t;
    while fn_0 < n_filters {
        *(*lower_freqs).data.offset(fn_0 as isize) =
            *(*freqs).data.offset(fn_0 as isize);
        *(*center_freqs).data.offset(fn_0 as isize) =
            *(*freqs).data.offset(fn_0.wrapping_add(1 as libc::c_int as
                                                        libc::c_uint) as
                                      isize);
        *(*upper_freqs).data.offset(fn_0 as isize) =
            *(*freqs).data.offset(fn_0.wrapping_add(2 as libc::c_int as
                                                        libc::c_uint) as
                                      isize);
        fn_0 = fn_0.wrapping_add(1)
    }
    /* compute triangle heights so that each triangle has unit area */
    if aubio_filterbank_get_norm(fb) != 0. {
        fn_0 = 0 as libc::c_int as uint_t;
        while fn_0 < n_filters {
            *(*triangle_heights).data.offset(fn_0 as isize) =
                (2.0f64 /
                     (*(*upper_freqs).data.offset(fn_0 as isize) -
                          *(*lower_freqs).data.offset(fn_0 as isize)) as
                         libc::c_double) as smpl_t;
            fn_0 = fn_0.wrapping_add(1)
        }
    } else { fvec_ones(triangle_heights); }
    /* fill fft_freqs lookup table, which assigns the frequency in hz to each bin */
    bin = 0 as libc::c_int as uint_t;
    while bin < win_s {
        *(*fft_freqs).data.offset(bin as isize) =
            aubio_bintofreq(bin as smpl_t, samplerate,
                            win_s.wrapping_sub(1 as libc::c_int as
                                                   libc::c_uint).wrapping_mul(2
                                                                                  as
                                                                                  libc::c_int
                                                                                  as
                                                                                  libc::c_uint)
                                as smpl_t);
        bin = bin.wrapping_add(1)
    }
    /* zeroing of all filters */
    fmat_zeros(filters);
    /* building each filter table */
    fn_0 = 0 as libc::c_int as uint_t;
    while fn_0 < n_filters {
        /* skip first elements */
        bin = 0 as libc::c_int as uint_t;
        while bin < win_s.wrapping_sub(1 as libc::c_int as libc::c_uint) {
            if *(*fft_freqs).data.offset(bin as isize) <=
                   *(*lower_freqs).data.offset(fn_0 as isize) &&
                   *(*fft_freqs).data.offset(bin.wrapping_add(1 as libc::c_int
                                                                  as
                                                                  libc::c_uint)
                                                 as isize) >
                       *(*lower_freqs).data.offset(fn_0 as isize) {
                bin = bin.wrapping_add(1);
                break ;
            } else { bin = bin.wrapping_add(1) }
        }
        /* nothing else to do */
        riseInc =
            *(*triangle_heights).data.offset(fn_0 as isize) /
                (*(*center_freqs).data.offset(fn_0 as isize) -
                     *(*lower_freqs).data.offset(fn_0 as isize));
        while bin < win_s.wrapping_sub(1 as libc::c_int as libc::c_uint)
              /* compute positive slope step size */
              /* compute coefficients in positive slope */
              {
            *(*(*filters).data.offset(fn_0 as isize)).offset(bin as isize) =
                (*(*fft_freqs).data.offset(bin as isize) -
                     *(*lower_freqs).data.offset(fn_0 as isize)) * riseInc;
            if *(*fft_freqs).data.offset(bin.wrapping_add(1 as libc::c_int as
                                                              libc::c_uint) as
                                             isize) >=
                   *(*center_freqs).data.offset(fn_0 as isize) {
                bin = bin.wrapping_add(1);
                break ;
            } else { bin = bin.wrapping_add(1) }
        }
        downInc =
            *(*triangle_heights).data.offset(fn_0 as isize) /
                (*(*upper_freqs).data.offset(fn_0 as isize) -
                     *(*center_freqs).data.offset(fn_0 as isize));
        while bin < win_s.wrapping_sub(1 as libc::c_int as libc::c_uint)
              /* compute negative slope step size */
              /* compute coefficents in negative slope */
              {
            let ref mut fresh0 =
                *(*(*filters).data.offset(fn_0 as
                                              isize)).offset(bin as isize);
            *fresh0 +=
                (*(*upper_freqs).data.offset(fn_0 as isize) -
                     *(*fft_freqs).data.offset(bin as isize)) * downInc;
            if (*(*(*filters).data.offset(fn_0 as isize)).offset(bin as isize)
                    as libc::c_double) < 0.0f64 {
                *(*(*filters).data.offset(fn_0 as isize)).offset(bin as isize)
                    = 0.0f64 as smpl_t
            }
            if *(*fft_freqs).data.offset(bin.wrapping_add(1 as libc::c_int as
                                                              libc::c_uint) as
                                             isize) >=
                   *(*upper_freqs).data.offset(fn_0 as isize) {
                break ;
            }
            bin = bin.wrapping_add(1)
        }
        fn_0 = fn_0.wrapping_add(1)
    }
    /* destroy temporarly allocated vectors */
    del_fvec(lower_freqs);
    del_fvec(upper_freqs);
    del_fvec(center_freqs);
    del_fvec(triangle_heights);
    del_fvec(fft_freqs);
    return AUBIO_OK as libc::c_int as uint_t;
}
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
pub unsafe extern "C" fn aubio_filterbank_set_mel_coeffs_slaney(mut fb:
                                                                    *mut aubio_filterbank_t,
                                                                mut samplerate:
                                                                    smpl_t)
 -> uint_t {
    /* Malcolm Slaney parameters */
    let lowestFrequency: smpl_t = 133.3333f64 as smpl_t;
    let linearSpacing: smpl_t = 66.66666666f64 as smpl_t;
    let logSpacing: smpl_t = 1.0711703f64 as smpl_t;
    let linearFilters: uint_t = 13 as libc::c_int as uint_t;
    let logFilters: uint_t = 27 as libc::c_int as uint_t;
    let n_filters: uint_t = linearFilters.wrapping_add(logFilters);
    let mut fn_0: uint_t = 0;
    let mut retval: uint_t = 0;
    let mut lastlinearCF: smpl_t = 0.;
    /* buffers to compute filter frequencies */
    let mut freqs: *mut fvec_t = 0 as *mut fvec_t;
    if samplerate <= 0 as libc::c_int as libc::c_float {
        aubio_log(AUBIO_LOG_ERR as libc::c_int,
                  b"AUBIO ERROR: filterbank: set_mel_coeffs_slaney samplerate should be > 0\n\x00"
                      as *const u8 as *const libc::c_char);
        return AUBIO_FAIL as libc::c_int as uint_t
    }
    freqs =
        new_fvec(n_filters.wrapping_add(2 as libc::c_int as libc::c_uint));
    /* first step: fill all the linear filter frequencies */
    fn_0 = 0 as libc::c_int as uint_t;
    while fn_0 < linearFilters {
        *(*freqs).data.offset(fn_0 as isize) =
            lowestFrequency + fn_0 as libc::c_float * linearSpacing;
        fn_0 = fn_0.wrapping_add(1)
    }
    lastlinearCF =
        *(*freqs).data.offset(fn_0.wrapping_sub(1 as libc::c_int as
                                                    libc::c_uint) as isize);
    /* second step: fill all the log filter frequencies */
    fn_0 = 0 as libc::c_int as uint_t;
    while fn_0 < logFilters.wrapping_add(2 as libc::c_int as libc::c_uint) {
        *(*freqs).data.offset(fn_0.wrapping_add(linearFilters) as isize) =
            lastlinearCF *
                powf(logSpacing,
                     fn_0.wrapping_add(1 as libc::c_int as libc::c_uint) as
                         libc::c_float);
        fn_0 = fn_0.wrapping_add(1)
    }
    /* now compute the actual coefficients */
    retval = aubio_filterbank_set_triangle_bands(fb, freqs, samplerate);
    /* destroy vector used to store frequency limits */
    del_fvec(freqs);
    return retval;
}
unsafe extern "C" fn aubio_filterbank_check_freqs(mut fb:
                                                      *mut aubio_filterbank_t,
                                                  mut samplerate: smpl_t,
                                                  mut freq_min: *mut smpl_t,
                                                  mut freq_max: *mut smpl_t)
 -> uint_t {
    if samplerate <= 0 as libc::c_int as libc::c_float {
        aubio_log(AUBIO_LOG_ERR as libc::c_int,
                  b"AUBIO ERROR: filterbank: set_mel_coeffs samplerate should be > 0\n\x00"
                      as *const u8 as *const libc::c_char);
        return AUBIO_FAIL as libc::c_int as uint_t
    }
    if *freq_max < 0 as libc::c_int as libc::c_float {
        aubio_log(AUBIO_LOG_ERR as libc::c_int,
                  b"AUBIO ERROR: filterbank: set_mel_coeffs freq_max should be > 0\n\x00"
                      as *const u8 as *const libc::c_char);
        return AUBIO_FAIL as libc::c_int as uint_t
    } else {
        if *freq_max == 0 as libc::c_int as libc::c_float {
            *freq_max = (samplerate as libc::c_double / 2.0f64) as smpl_t
        }
    }
    if *freq_min < 0 as libc::c_int as libc::c_float {
        aubio_log(AUBIO_LOG_ERR as libc::c_int,
                  b"AUBIO ERROR: filterbank: set_mel_coeffs freq_min should be > 0\n\x00"
                      as *const u8 as *const libc::c_char);
        return AUBIO_FAIL as libc::c_int as uint_t
    }
    return AUBIO_OK as libc::c_int as uint_t;
}
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
pub unsafe extern "C" fn aubio_filterbank_set_mel_coeffs(mut fb:
                                                             *mut aubio_filterbank_t,
                                                         mut samplerate:
                                                             smpl_t,
                                                         mut freq_min: smpl_t,
                                                         mut freq_max: smpl_t)
 -> uint_t {
    let mut m: uint_t = 0;
    let mut retval: uint_t = 0;
    let mut start: smpl_t = freq_min;
    let mut end: smpl_t = freq_max;
    let mut step: smpl_t = 0.;
    let mut freqs: *mut fvec_t = 0 as *mut fvec_t;
    let mut coeffs: *mut fmat_t = aubio_filterbank_get_coeffs(fb);
    let mut n_bands: uint_t = (*coeffs).height;
    if aubio_filterbank_check_freqs(fb, samplerate, &mut start, &mut end) != 0
       {
        return AUBIO_FAIL as libc::c_int as uint_t
    }
    start = aubio_hztomel(start);
    end = aubio_hztomel(end);
    freqs = new_fvec(n_bands.wrapping_add(2 as libc::c_int as libc::c_uint));
    step =
        (end - start) /
            n_bands.wrapping_add(1 as libc::c_int as libc::c_uint) as
                libc::c_float;
    m = 0 as libc::c_int as uint_t;
    while m < n_bands.wrapping_add(2 as libc::c_int as libc::c_uint) {
        *(*freqs).data.offset(m as isize) =
            if (aubio_meltohz(start + step * m as libc::c_float) as
                    libc::c_double) < samplerate as libc::c_double / 2.0f64 {
                aubio_meltohz(start + step * m as libc::c_float) as
                    libc::c_double
            } else { (samplerate as libc::c_double) / 2.0f64 } as smpl_t;
        m = m.wrapping_add(1)
    }
    retval = aubio_filterbank_set_triangle_bands(fb, freqs, samplerate);
    /* destroy vector used to store frequency limits */
    del_fvec(freqs);
    return retval;
}
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
pub unsafe extern "C" fn aubio_filterbank_set_mel_coeffs_htk(mut fb:
                                                                 *mut aubio_filterbank_t,
                                                             mut samplerate:
                                                                 smpl_t,
                                                             mut freq_min:
                                                                 smpl_t,
                                                             mut freq_max:
                                                                 smpl_t)
 -> uint_t {
    let mut m: uint_t = 0;
    let mut retval: uint_t = 0;
    let mut start: smpl_t = freq_min;
    let mut end: smpl_t = freq_max;
    let mut step: smpl_t = 0.;
    let mut freqs: *mut fvec_t = 0 as *mut fvec_t;
    let mut coeffs: *mut fmat_t = aubio_filterbank_get_coeffs(fb);
    let mut n_bands: uint_t = (*coeffs).height;
    if aubio_filterbank_check_freqs(fb, samplerate, &mut start, &mut end) != 0
       {
        return AUBIO_FAIL as libc::c_int as uint_t
    }
    start = aubio_hztomel_htk(start);
    end = aubio_hztomel_htk(end);
    freqs = new_fvec(n_bands.wrapping_add(2 as libc::c_int as libc::c_uint));
    step =
        (end - start) /
            n_bands.wrapping_add(1 as libc::c_int as libc::c_uint) as
                libc::c_float;
    m = 0 as libc::c_int as uint_t;
    while m < n_bands.wrapping_add(2 as libc::c_int as libc::c_uint) {
        *(*freqs).data.offset(m as isize) =
            if (aubio_meltohz_htk(start + step * m as libc::c_float) as
                    libc::c_double) < samplerate as libc::c_double / 2.0f64 {
                aubio_meltohz_htk(start + step * m as libc::c_float) as
                    libc::c_double
            } else { (samplerate as libc::c_double) / 2.0f64 } as smpl_t;
        m = m.wrapping_add(1)
    }
    retval = aubio_filterbank_set_triangle_bands(fb, freqs, samplerate);
    /* destroy vector used to store frequency limits */
    del_fvec(freqs);
    return retval;
}
