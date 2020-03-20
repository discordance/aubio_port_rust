use crate::transpiled::hist::aubio_hist_t;

extern "C" {
    // pub type _aubio_hist_t;
    #[no_mangle]
    fn cosf(_: f32) -> f32;
    #[no_mangle]
    fn aubio_log(level: sint_t, fmt: *const char_t, _: ...) -> uint_t;
    #[no_mangle]
    fn calloc(_: u64, _: u64) -> *mut core::ffi::c_void;
    #[no_mangle]
    fn free(_: *mut core::ffi::c_void);
    #[no_mangle]
    fn logf(_: f32) -> f32;
    #[no_mangle]
    fn fabsf(_: f32) -> f32;
    #[no_mangle]
    fn sqrtf(_: f32) -> f32;
    #[no_mangle]
    fn strcmp(_: *const i8, _: *const i8) -> i32;
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
    /* * compute the principal argument

      This function maps the input phase to its corresponding value wrapped in the
    range \f$ [-\pi, \pi] \f$.

      \param phase unwrapped phase to map to the unit circle

      \return equivalent phase wrapped to the unit circle

    */
    #[no_mangle]
    fn aubio_unwrap2pi(phase: smpl_t) -> smpl_t;
    /* * histogram creation

      \param flow minimum input
      \param fhig maximum input
      \param nelems number of histogram columns

    */
    #[no_mangle]
    fn new_aubio_hist(flow: smpl_t, fhig: smpl_t, nelems: uint_t) -> *mut aubio_hist_t;
    /* * histogram deletion */
    #[no_mangle]
    fn del_aubio_hist(s: *mut aubio_hist_t);
    /* * compute the mean of the histogram */
    #[no_mangle]
    fn aubio_hist_mean(s: *const aubio_hist_t) -> smpl_t;
    /* * weight the histogram */
    #[no_mangle]
    fn aubio_hist_weight(s: *mut aubio_hist_t);
    /* * compute dynamic histogram for non-null elements */
    #[no_mangle]
    fn aubio_hist_dyn_notnull(s: *mut aubio_hist_t, input: *mut fvec_t);
    #[no_mangle]
    fn aubio_specdesc_centroid(o: *mut aubio_specdesc_t, spec: *const cvec_t, desc: *mut fvec_t);
    #[no_mangle]
    fn aubio_specdesc_spread(o: *mut aubio_specdesc_t, spec: *const cvec_t, desc: *mut fvec_t);
    #[no_mangle]
    fn aubio_specdesc_skewness(o: *mut aubio_specdesc_t, spec: *const cvec_t, desc: *mut fvec_t);
    #[no_mangle]
    fn aubio_specdesc_kurtosis(o: *mut aubio_specdesc_t, spec: *const cvec_t, desc: *mut fvec_t);
    #[no_mangle]
    fn aubio_specdesc_slope(o: *mut aubio_specdesc_t, spec: *const cvec_t, desc: *mut fvec_t);
    #[no_mangle]
    fn aubio_specdesc_decrease(o: *mut aubio_specdesc_t, spec: *const cvec_t, desc: *mut fvec_t);
    #[no_mangle]
    fn aubio_specdesc_rolloff(o: *mut aubio_specdesc_t, spec: *const cvec_t, desc: *mut fvec_t);
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
pub type smpl_t = f32;
/* * print format for sample in double precision */
/* * unsigned integer */
pub type uint_t = u32;
/* * signed integer */
pub type sint_t = i32;
/* * character */
pub type char_t = i8;
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
pub type aubio_log_level = u32;
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
/* * structure to store object state */
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _aubio_specdesc_t {
    pub onset_type: aubio_specdesc_type,
    pub funcpointer: Option<
        unsafe extern "C" fn(_: *mut aubio_specdesc_t, _: *const cvec_t, _: *mut fvec_t) -> (),
    >,
    pub threshold: smpl_t,
    pub oldmag: *mut fvec_t,
    pub dev1: *mut fvec_t,
    pub theta1: *mut fvec_t,
    pub theta2: *mut fvec_t,
    pub histog: *mut aubio_hist_t,
}
/*
  Copyright (C) 2003-2013 Paul Brossier <piem@aubio.org>

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
/*
  Copyright (C) 2003-2013 Paul Brossier <piem@aubio.org>

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

  Spectral description functions

  All of the following spectral description functions take as arguments the FFT
  of a windowed signal (as created with aubio_pvoc). They output one smpl_t per
  buffer (stored in a vector of size [1]).

  \section specdesc Spectral description functions

  A list of the spectral description methods currently available follows.

  \subsection onsetdesc Onset detection functions

  These functions are designed to raise at notes attacks in music signals.

  \b \p energy : Energy based onset detection function

  This function calculates the local energy of the input spectral frame.

  \b \p hfc : High Frequency Content onset detection function

  This method computes the High Frequency Content (HFC) of the input spectral
  frame. The resulting function is efficient at detecting percussive onsets.

  Paul Masri. Computer modeling of Sound for Transformation and Synthesis of
  Musical Signal. PhD dissertation, University of Bristol, UK, 1996.

  \b \p complex : Complex Domain Method onset detection function

  Christopher Duxbury, Mike E. Davies, and Mark B. Sandler. Complex domain
  onset detection for musical signals. In Proceedings of the Digital Audio
  Effects Conference, DAFx-03, pages 90-93, London, UK, 2003.

  \b \p phase : Phase Based Method onset detection function

  Juan-Pablo Bello, Mike P. Davies, and Mark B. Sandler. Phase-based note onset
  detection for music signals. In Proceedings of the IEEE International
  Conference on Acoustics Speech and Signal Processing, pages 441­444,
  Hong-Kong, 2003.

  \b \p wphase : Weighted Phase Deviation onset detection function

  S. Dixon. Onset detection revisited. In Proceedings of the 9th International
  Conference on Digital Audio Ef- fects (DAFx) , pages 133–137, 2006.

  http://www.eecs.qmul.ac.uk/~simond/pub/2006/dafx.pdf

  \b \p specdiff : Spectral difference method onset detection function

  Jonhatan Foote and Shingo Uchihashi. The beat spectrum: a new approach to
  rhythm analysis. In IEEE International Conference on Multimedia and Expo
  (ICME 2001), pages 881­884, Tokyo, Japan, August 2001.

  \b \p kl : Kullback-Liebler onset detection function

  Stephen Hainsworth and Malcom Macleod. Onset detection in music audio
  signals. In Proceedings of the International Computer Music Conference
  (ICMC), Singapore, 2003.

  \b \p mkl : Modified Kullback-Liebler onset detection function

  Paul Brossier, ``Automatic annotation of musical audio for interactive
  systems'', Chapter 2, Temporal segmentation, PhD thesis, Centre for Digital
  music, Queen Mary University of London, London, UK, 2006.

  \b \p specflux : Spectral Flux

  Simon Dixon, Onset Detection Revisited, in ``Proceedings of the 9th
  International Conference on Digital Audio Effects'' (DAFx-06), Montreal,
  Canada, 2006.

  \subsection shapedesc Spectral shape descriptors

  The following descriptors are described in:

  Geoffroy Peeters, <i>A large set of audio features for sound description
  (similarity and classification) in the CUIDADO project</i>, CUIDADO I.S.T.
  Project Report 2004 (<a
  href="http://www.ircam.fr/anasyn/peeters/ARTICLES/Peeters_2003_cuidadoaudiofeatures.pdf">pdf</a>)

  \b \p centroid : Spectral centroid

  The spectral centroid represents the barycenter of the spectrum.

  \e Note: This function returns the result in bin. To get the spectral
  centroid in Hz, aubio_bintofreq() should be used.

  \b \p spread : Spectral spread

  The spectral spread is the variance of the spectral distribution around its
  centroid.

  See also <a href="http://en.wikipedia.org/wiki/Standard_deviation">Standard
  deviation</a> on Wikipedia.

  \b \p skewness : Spectral skewness

  Similarly, the skewness is computed from the third order moment of the
  spectrum. A negative skewness indicates more energy on the lower part of the
  spectrum. A positive skewness indicates more energy on the high frequency of
  the spectrum.

  See also <a href="http://en.wikipedia.org/wiki/Skewness">Skewness</a> on
  Wikipedia.

  \b \p kurtosis : Spectral kurtosis

  The kurtosis is a measure of the flatness of the spectrum, computed from the
  fourth order moment.

  See also <a href="http://en.wikipedia.org/wiki/Kurtosis">Kurtosis</a> on
  Wikipedia.

  \b \p slope : Spectral slope

  The spectral slope represents decreasing rate of the spectral amplitude,
  computed using a linear regression.

  \b \p decrease : Spectral decrease

  The spectral decrease is another representation of the decreasing rate,
  based on perceptual criteria.

  \b \p rolloff : Spectral roll-off

  This function returns the bin number below which 95% of the spectrum energy
  is found.

  \example spectral/test-specdesc.c

*/
/* * spectral description structure */
pub type aubio_specdesc_t = _aubio_specdesc_t;
pub type aubio_specdesc_type = u32;
pub const aubio_onset_default: aubio_specdesc_type = 2;
pub const aubio_specmethod_rolloff: aubio_specdesc_type = 15;
pub const aubio_specmethod_decrease: aubio_specdesc_type = 14;
pub const aubio_specmethod_slope: aubio_specdesc_type = 13;
pub const aubio_specmethod_kurtosis: aubio_specdesc_type = 12;
pub const aubio_specmethod_skewness: aubio_specdesc_type = 11;
pub const aubio_specmethod_spread: aubio_specdesc_type = 10;
pub const aubio_specmethod_centroid: aubio_specdesc_type = 9;
pub const aubio_onset_specflux: aubio_specdesc_type = 8;
pub const aubio_onset_mkl: aubio_specdesc_type = 7;
pub const aubio_onset_kl: aubio_specdesc_type = 6;
pub const aubio_onset_wphase: aubio_specdesc_type = 5;
pub const aubio_onset_phase: aubio_specdesc_type = 4;
pub const aubio_onset_complex: aubio_specdesc_type = 3;
pub const aubio_onset_hfc: aubio_specdesc_type = 2;
pub const aubio_onset_specdiff: aubio_specdesc_type = 1;
pub const aubio_onset_energy: aubio_specdesc_type = 0;
#[inline(always)]
unsafe extern "C" fn __inline_isnanf(mut __x: f32) -> i32 {
    return (__x != __x) as i32;
}
#[inline(always)]
unsafe extern "C" fn __inline_isnand(mut __x: f64) -> i32 {
    return (__x != __x) as i32;
}
#[inline(always)]
// rust does not support f128
unsafe extern "C" fn __inline_isnanl(mut __x: f64) -> i32 {
    return (__x != __x) as i32;
}
/*
  Copyright (C) 2003-2009 Paul Brossier <piem@aubio.org>

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
/* *< histogram */
/* Energy based onset detection function */
#[no_mangle]
pub unsafe extern "C" fn aubio_specdesc_energy(
    mut o: *mut aubio_specdesc_t,
    mut fftgrain: *const cvec_t,
    mut onset: *mut fvec_t,
) {
    let mut j: uint_t = 0;
    *(*onset).data.offset(0 as i32 as isize) = 0.0f64 as smpl_t;
    j = 0 as i32 as uint_t;
    while j < (*fftgrain).length {
        let ref mut fresh0 = *(*onset).data.offset(0 as i32 as isize);
        *fresh0 += *(*fftgrain).norm.offset(j as isize) * *(*fftgrain).norm.offset(j as isize);
        j = j.wrapping_add(1)
    }
}
/* High Frequency Content onset detection function */
#[no_mangle]
pub unsafe extern "C" fn aubio_specdesc_hfc(
    mut o: *mut aubio_specdesc_t,
    mut fftgrain: *const cvec_t,
    mut onset: *mut fvec_t,
) {
    let mut j: uint_t = 0;
    *(*onset).data.offset(0 as i32 as isize) = 0.0f64 as smpl_t;
    j = 0 as i32 as uint_t;
    while j < (*fftgrain).length {
        let ref mut fresh1 = *(*onset).data.offset(0 as i32 as isize);
        *fresh1 += j.wrapping_add(1 as i32 as u32) as f32 * *(*fftgrain).norm.offset(j as isize);
        j = j.wrapping_add(1)
    }
}
/* Complex Domain Method onset detection function */
#[no_mangle]
pub unsafe extern "C" fn aubio_specdesc_complex(
    mut o: *mut aubio_specdesc_t,
    mut fftgrain: *const cvec_t,
    mut onset: *mut fvec_t,
) {
    let mut j: uint_t = 0;
    let mut nbins: uint_t = (*fftgrain).length;
    *(*onset).data.offset(0 as i32 as isize) = 0.0f64 as smpl_t;
    j = 0 as i32 as uint_t;
    while j < nbins {
        // compute the predicted phase
        *(*(*o).dev1).data.offset(j as isize) =
            (2.0f64 * *(*(*o).theta1).data.offset(j as isize) as f64
                - *(*(*o).theta2).data.offset(j as isize) as f64) as smpl_t;
        // compute the euclidean distance in the complex domain
        // sqrt ( r_1^2 + r_2^2 - 2 * r_1 * r_2 * \cos ( \phi_1 - \phi_2 ) )
        let ref mut fresh2 = *(*onset).data.offset(0 as i32 as isize);
        *fresh2 += sqrtf(fabsf(
            *(*(*o).oldmag).data.offset(j as isize) * *(*(*o).oldmag).data.offset(j as isize)
                + *(*fftgrain).norm.offset(j as isize) * *(*fftgrain).norm.offset(j as isize)
                - 2 as i32 as f32
                    * *(*(*o).oldmag).data.offset(j as isize)
                    * *(*fftgrain).norm.offset(j as isize)
                    * cosf(
                        *(*(*o).dev1).data.offset(j as isize)
                            - *(*fftgrain).phas.offset(j as isize),
                    ),
        ));
        /* swap old phase data (need to remember 2 frames behind)*/
        *(*(*o).theta2).data.offset(j as isize) = *(*(*o).theta1).data.offset(j as isize);
        *(*(*o).theta1).data.offset(j as isize) = *(*fftgrain).phas.offset(j as isize);
        /* swap old magnitude data (1 frame is enough) */
        *(*(*o).oldmag).data.offset(j as isize) = *(*fftgrain).norm.offset(j as isize);
        j = j.wrapping_add(1)
    }
}
/* Phase Based Method onset detection function */
#[no_mangle]
pub unsafe extern "C" fn aubio_specdesc_phase(
    mut o: *mut aubio_specdesc_t,
    mut fftgrain: *const cvec_t,
    mut onset: *mut fvec_t,
) {
    let mut j: uint_t = 0;
    let mut nbins: uint_t = (*fftgrain).length;
    *(*onset).data.offset(0 as i32 as isize) = 0.0f64 as smpl_t;
    *(*(*o).dev1).data.offset(0 as i32 as isize) = 0.0f64 as smpl_t;
    j = 0 as i32 as uint_t;
    while j < nbins {
        *(*(*o).dev1).data.offset(j as isize) = aubio_unwrap2pi(
            (*(*fftgrain).phas.offset(j as isize) as f64
                - 2.0f64 * *(*(*o).theta1).data.offset(j as isize) as f64
                + *(*(*o).theta2).data.offset(j as isize) as f64) as smpl_t,
        );
        if (*o).threshold < *(*fftgrain).norm.offset(j as isize) {
            *(*(*o).dev1).data.offset(j as isize) = fabsf(*(*(*o).dev1).data.offset(j as isize))
        } else {
            *(*(*o).dev1).data.offset(j as isize) = 0.0f64 as smpl_t
        }
        /* keep a track of the past frames */
        *(*(*o).theta2).data.offset(j as isize) = *(*(*o).theta1).data.offset(j as isize);
        *(*(*o).theta1).data.offset(j as isize) = *(*fftgrain).phas.offset(j as isize);
        j = j.wrapping_add(1)
    }
    /* apply o->histogram */
    aubio_hist_dyn_notnull((*o).histog, (*o).dev1);
    /* weight it */
    aubio_hist_weight((*o).histog);
    /* its mean is the result */
    *(*onset).data.offset(0 as i32 as isize) = aubio_hist_mean((*o).histog);
    //onset->data[0] = fvec_mean(o->dev1);
}
/* weighted phase */
#[no_mangle]
pub unsafe extern "C" fn aubio_specdesc_wphase(
    mut o: *mut aubio_specdesc_t,
    mut fftgrain: *const cvec_t,
    mut onset: *mut fvec_t,
) {
    let mut i: uint_t = 0;
    aubio_specdesc_phase(o, fftgrain, onset);
    i = 0 as i32 as uint_t;
    while i < (*fftgrain).length {
        let ref mut fresh3 = *(*(*o).dev1).data.offset(i as isize);
        *fresh3 *= *(*fftgrain).norm.offset(i as isize);
        i = i.wrapping_add(1)
    }
    /* apply o->histogram */
    aubio_hist_dyn_notnull((*o).histog, (*o).dev1);
    /* weight it */
    aubio_hist_weight((*o).histog);
    /* its mean is the result */
    *(*onset).data.offset(0 as i32 as isize) = aubio_hist_mean((*o).histog);
}
/* Spectral difference method onset detection function */
#[no_mangle]
pub unsafe extern "C" fn aubio_specdesc_specdiff(
    mut o: *mut aubio_specdesc_t,
    mut fftgrain: *const cvec_t,
    mut onset: *mut fvec_t,
) {
    let mut j: uint_t = 0;
    let mut nbins: uint_t = (*fftgrain).length;
    *(*onset).data.offset(0 as i32 as isize) = 0.0f64 as smpl_t;
    j = 0 as i32 as uint_t;
    while j < nbins {
        *(*(*o).dev1).data.offset(j as isize) = sqrtf(fabsf(
            *(*fftgrain).norm.offset(j as isize) * *(*fftgrain).norm.offset(j as isize)
                - *(*(*o).oldmag).data.offset(j as isize) * *(*(*o).oldmag).data.offset(j as isize),
        ));
        if (*o).threshold < *(*fftgrain).norm.offset(j as isize) {
            *(*(*o).dev1).data.offset(j as isize) = fabsf(*(*(*o).dev1).data.offset(j as isize))
        } else {
            *(*(*o).dev1).data.offset(j as isize) = 0.0f64 as smpl_t
        }
        *(*(*o).oldmag).data.offset(j as isize) = *(*fftgrain).norm.offset(j as isize);
        j = j.wrapping_add(1)
    }
    /* apply o->histogram (act somewhat as a low pass on the
     * overall function)*/
    aubio_hist_dyn_notnull((*o).histog, (*o).dev1);
    /* weight it */
    aubio_hist_weight((*o).histog);
    /* its mean is the result */
    *(*onset).data.offset(0 as i32 as isize) = aubio_hist_mean((*o).histog);
}
/* Kullback Liebler onset detection function
 * note we use ln(1+Xn/(Xn-1+0.0001)) to avoid
 * negative (1.+) and infinite values (+1.e-10) */
#[no_mangle]
pub unsafe extern "C" fn aubio_specdesc_kl(
    mut o: *mut aubio_specdesc_t,
    mut fftgrain: *const cvec_t,
    mut onset: *mut fvec_t,
) {
    let mut j: uint_t = 0;
    *(*onset).data.offset(0 as i32 as isize) = 0.0f64 as smpl_t;
    j = 0 as i32 as uint_t;
    while j < (*fftgrain).length {
        let ref mut fresh4 = *(*onset).data.offset(0 as i32 as isize);
        *fresh4 += *(*fftgrain).norm.offset(j as isize)
            * logf(
                (1.0f64
                    + *(*fftgrain).norm.offset(j as isize) as f64
                        / (*(*(*o).oldmag).data.offset(j as isize) as f64 + 1.0e-1f64))
                    as f32,
            );
        *(*(*o).oldmag).data.offset(j as isize) = *(*fftgrain).norm.offset(j as isize);
        j = j.wrapping_add(1)
    }
    if if ::std::mem::size_of::<smpl_t>() as u64 == ::std::mem::size_of::<f32>() as u64 {
        __inline_isnanf(*(*onset).data.offset(0 as i32 as isize))
    } else if ::std::mem::size_of::<smpl_t>() as u64 == ::std::mem::size_of::<f64>() as u64 {
        __inline_isnand(*(*onset).data.offset(0 as i32 as isize) as f64)
    } else {
        __inline_isnand(*(*onset).data.offset(0 as i32 as isize) as f64)
    } != 0
    {
        *(*onset).data.offset(0 as i32 as isize) = 0.0f64 as smpl_t
    };
}
/* Modified Kullback Liebler onset detection function
 * note we use ln(1+Xn/(Xn-1+0.0001)) to avoid
 * negative (1.+) and infinite values (+1.e-10) */
#[no_mangle]
pub unsafe extern "C" fn aubio_specdesc_mkl(
    mut o: *mut aubio_specdesc_t,
    mut fftgrain: *const cvec_t,
    mut onset: *mut fvec_t,
) {
    let mut j: uint_t = 0;
    *(*onset).data.offset(0 as i32 as isize) = 0.0f64 as smpl_t;
    j = 0 as i32 as uint_t;
    while j < (*fftgrain).length {
        let ref mut fresh5 = *(*onset).data.offset(0 as i32 as isize);
        *fresh5 += logf(
            (1.0f64
                + *(*fftgrain).norm.offset(j as isize) as f64
                    / (*(*(*o).oldmag).data.offset(j as isize) as f64 + 1.0e-1f64))
                as f32,
        );
        *(*(*o).oldmag).data.offset(j as isize) = *(*fftgrain).norm.offset(j as isize);
        j = j.wrapping_add(1)
    }
    if if ::std::mem::size_of::<smpl_t>() as u64 == ::std::mem::size_of::<f32>() as u64 {
        __inline_isnanf(*(*onset).data.offset(0 as i32 as isize))
    } else if ::std::mem::size_of::<smpl_t>() as u64 == ::std::mem::size_of::<f64>() as u64 {
        __inline_isnand(*(*onset).data.offset(0 as i32 as isize) as f64)
    } else {
        __inline_isnand(*(*onset).data.offset(0 as i32 as isize) as f64)
    } != 0
    {
        *(*onset).data.offset(0 as i32 as isize) = 0.0f64 as smpl_t
    };
}
/* Spectral flux */
#[no_mangle]
pub unsafe extern "C" fn aubio_specdesc_specflux(
    mut o: *mut aubio_specdesc_t,
    mut fftgrain: *const cvec_t,
    mut onset: *mut fvec_t,
) {
    let mut j: uint_t = 0;
    *(*onset).data.offset(0 as i32 as isize) = 0.0f64 as smpl_t;
    j = 0 as i32 as uint_t;
    while j < (*fftgrain).length {
        if *(*fftgrain).norm.offset(j as isize) > *(*(*o).oldmag).data.offset(j as isize) {
            let ref mut fresh6 = *(*onset).data.offset(0 as i32 as isize);
            *fresh6 +=
                *(*fftgrain).norm.offset(j as isize) - *(*(*o).oldmag).data.offset(j as isize)
        }
        *(*(*o).oldmag).data.offset(j as isize) = *(*fftgrain).norm.offset(j as isize);
        j = j.wrapping_add(1)
    }
}
/* * execute spectral description function on a spectral frame

  Generic function to compute spectral description.

  \param o spectral description object as returned by new_aubio_specdesc()
  \param fftgrain input signal spectrum as computed by aubio_pvoc_do
  \param desc output vector (one sample long, to send to the peak picking)

*/
/* Generic function pointing to the choosen one */
#[no_mangle]
pub unsafe extern "C" fn aubio_specdesc_do(
    mut o: *mut aubio_specdesc_t,
    mut fftgrain: *const cvec_t,
    mut onset: *mut fvec_t,
) {
    (*o).funcpointer.expect("non-null function pointer")(o, fftgrain, onset);
}
/* * creation of a spectral description object

  \param method spectral description method
  \param buf_size length of the input spectrum frame

  The parameter \p method is a string that can be any of:

    - onset novelty functions: `complex`, `energy`, `hfc`, `kl`, `mkl`,
    `phase`, `specdiff`, `specflux`, `wphase`,

    - spectral descriptors: `centroid`, `decrease`, `kurtosis`, `rolloff`,
    `skewness`, `slope`, `spread`.

*/
/* Allocate memory for an onset detection
 * depending on the choosen type, allocate memory as needed
 */
#[no_mangle]
pub unsafe extern "C" fn new_aubio_specdesc(
    mut onset_mode: *const char_t,
    mut size: uint_t,
) -> *mut aubio_specdesc_t {
    let mut o: *mut aubio_specdesc_t = calloc(
        ::std::mem::size_of::<aubio_specdesc_t>() as u64,
        1 as i32 as u64,
    ) as *mut aubio_specdesc_t;
    let mut rsize: uint_t = size
        .wrapping_div(2 as i32 as u32)
        .wrapping_add(1 as i32 as u32);
    let mut onset_type: aubio_specdesc_type = aubio_onset_energy;
    if strcmp(onset_mode, b"energy\x00" as *const u8 as *const i8) == 0 as i32 {
        onset_type = aubio_onset_energy
    } else if strcmp(onset_mode, b"specdiff\x00" as *const u8 as *const i8) == 0 as i32 {
        onset_type = aubio_onset_specdiff
    } else if strcmp(onset_mode, b"hfc\x00" as *const u8 as *const i8) == 0 as i32 {
        onset_type = aubio_onset_hfc
    } else if strcmp(onset_mode, b"complexdomain\x00" as *const u8 as *const i8) == 0 as i32 {
        onset_type = aubio_onset_complex
    } else if strcmp(onset_mode, b"complex\x00" as *const u8 as *const i8) == 0 as i32 {
        onset_type = aubio_onset_complex
    } else if strcmp(onset_mode, b"phase\x00" as *const u8 as *const i8) == 0 as i32 {
        onset_type = aubio_onset_phase
    } else if strcmp(onset_mode, b"wphase\x00" as *const u8 as *const i8) == 0 as i32 {
        onset_type = aubio_onset_wphase
    } else if strcmp(onset_mode, b"mkl\x00" as *const u8 as *const i8) == 0 as i32 {
        onset_type = aubio_onset_mkl
    } else if strcmp(onset_mode, b"kl\x00" as *const u8 as *const i8) == 0 as i32 {
        onset_type = aubio_onset_kl
    } else if strcmp(onset_mode, b"specflux\x00" as *const u8 as *const i8) == 0 as i32 {
        onset_type = aubio_onset_specflux
    } else if strcmp(onset_mode, b"centroid\x00" as *const u8 as *const i8) == 0 as i32 {
        onset_type = aubio_specmethod_centroid
    } else if strcmp(onset_mode, b"spread\x00" as *const u8 as *const i8) == 0 as i32 {
        onset_type = aubio_specmethod_spread
    } else if strcmp(onset_mode, b"skewness\x00" as *const u8 as *const i8) == 0 as i32 {
        onset_type = aubio_specmethod_skewness
    } else if strcmp(onset_mode, b"kurtosis\x00" as *const u8 as *const i8) == 0 as i32 {
        onset_type = aubio_specmethod_kurtosis
    } else if strcmp(onset_mode, b"slope\x00" as *const u8 as *const i8) == 0 as i32 {
        onset_type = aubio_specmethod_slope
    } else if strcmp(onset_mode, b"decrease\x00" as *const u8 as *const i8) == 0 as i32 {
        onset_type = aubio_specmethod_decrease
    } else if strcmp(onset_mode, b"rolloff\x00" as *const u8 as *const i8) == 0 as i32 {
        onset_type = aubio_specmethod_rolloff
    } else if strcmp(onset_mode, b"old_default\x00" as *const u8 as *const i8) == 0 as i32 {
        onset_type = aubio_onset_default
    } else if strcmp(onset_mode, b"default\x00" as *const u8 as *const i8) == 0 as i32 {
        onset_type = aubio_onset_default
    } else {
        aubio_log(
            AUBIO_LOG_ERR as i32,
            b"AUBIO ERROR: specdesc: unknown spectral descriptor type \'%s\'\n\x00" as *const u8
                as *const i8,
            onset_mode,
        );
        free(o as *mut core::ffi::c_void);
        return 0 as *mut aubio_specdesc_t;
    }
    match onset_type as u32 {
        3 => {
            /* the other approaches will need some more memory spaces */
            (*o).oldmag = new_fvec(rsize);
            (*o).dev1 = new_fvec(rsize);
            (*o).theta1 = new_fvec(rsize);
            (*o).theta2 = new_fvec(rsize)
        }
        4 | 5 => {
            (*o).dev1 = new_fvec(rsize);
            (*o).theta1 = new_fvec(rsize);
            (*o).theta2 = new_fvec(rsize);
            (*o).histog = new_aubio_hist(
                0.0f64 as smpl_t,
                3.14159265358979323846264338327950288f64 as smpl_t,
                10 as i32 as uint_t,
            );
            (*o).threshold = 0.1f64 as smpl_t
        }
        1 => {
            (*o).oldmag = new_fvec(rsize);
            (*o).dev1 = new_fvec(rsize);
            (*o).histog = new_aubio_hist(
                0.0f64 as smpl_t,
                3.14159265358979323846264338327950288f64 as smpl_t,
                10 as i32 as uint_t,
            );
            (*o).threshold = 0.1f64 as smpl_t
        }
        6 | 7 | 8 => (*o).oldmag = new_fvec(rsize),
        0 | 2 | _ => {}
    }
    match onset_type as u32 {
        0 => {
            (*o).funcpointer = Some(
                aubio_specdesc_energy
                    as unsafe extern "C" fn(
                        _: *mut aubio_specdesc_t,
                        _: *const cvec_t,
                        _: *mut fvec_t,
                    ) -> (),
            )
        }
        2 => {
            (*o).funcpointer = Some(
                aubio_specdesc_hfc
                    as unsafe extern "C" fn(
                        _: *mut aubio_specdesc_t,
                        _: *const cvec_t,
                        _: *mut fvec_t,
                    ) -> (),
            )
        }
        3 => {
            (*o).funcpointer = Some(
                aubio_specdesc_complex
                    as unsafe extern "C" fn(
                        _: *mut aubio_specdesc_t,
                        _: *const cvec_t,
                        _: *mut fvec_t,
                    ) -> (),
            )
        }
        4 => {
            (*o).funcpointer = Some(
                aubio_specdesc_phase
                    as unsafe extern "C" fn(
                        _: *mut aubio_specdesc_t,
                        _: *const cvec_t,
                        _: *mut fvec_t,
                    ) -> (),
            )
        }
        5 => {
            (*o).funcpointer = Some(
                aubio_specdesc_wphase
                    as unsafe extern "C" fn(
                        _: *mut aubio_specdesc_t,
                        _: *const cvec_t,
                        _: *mut fvec_t,
                    ) -> (),
            )
        }
        1 => {
            (*o).funcpointer = Some(
                aubio_specdesc_specdiff
                    as unsafe extern "C" fn(
                        _: *mut aubio_specdesc_t,
                        _: *const cvec_t,
                        _: *mut fvec_t,
                    ) -> (),
            )
        }
        6 => {
            (*o).funcpointer = Some(
                aubio_specdesc_kl
                    as unsafe extern "C" fn(
                        _: *mut aubio_specdesc_t,
                        _: *const cvec_t,
                        _: *mut fvec_t,
                    ) -> (),
            )
        }
        7 => {
            (*o).funcpointer = Some(
                aubio_specdesc_mkl
                    as unsafe extern "C" fn(
                        _: *mut aubio_specdesc_t,
                        _: *const cvec_t,
                        _: *mut fvec_t,
                    ) -> (),
            )
        }
        8 => {
            (*o).funcpointer = Some(
                aubio_specdesc_specflux
                    as unsafe extern "C" fn(
                        _: *mut aubio_specdesc_t,
                        _: *const cvec_t,
                        _: *mut fvec_t,
                    ) -> (),
            )
        }
        9 => {
            (*o).funcpointer = Some(
                aubio_specdesc_centroid
                    as unsafe extern "C" fn(
                        _: *mut aubio_specdesc_t,
                        _: *const cvec_t,
                        _: *mut fvec_t,
                    ) -> (),
            )
        }
        10 => {
            (*o).funcpointer = Some(
                aubio_specdesc_spread
                    as unsafe extern "C" fn(
                        _: *mut aubio_specdesc_t,
                        _: *const cvec_t,
                        _: *mut fvec_t,
                    ) -> (),
            )
        }
        11 => {
            (*o).funcpointer = Some(
                aubio_specdesc_skewness
                    as unsafe extern "C" fn(
                        _: *mut aubio_specdesc_t,
                        _: *const cvec_t,
                        _: *mut fvec_t,
                    ) -> (),
            )
        }
        12 => {
            (*o).funcpointer = Some(
                aubio_specdesc_kurtosis
                    as unsafe extern "C" fn(
                        _: *mut aubio_specdesc_t,
                        _: *const cvec_t,
                        _: *mut fvec_t,
                    ) -> (),
            )
        }
        13 => {
            (*o).funcpointer = Some(
                aubio_specdesc_slope
                    as unsafe extern "C" fn(
                        _: *mut aubio_specdesc_t,
                        _: *const cvec_t,
                        _: *mut fvec_t,
                    ) -> (),
            )
        }
        14 => {
            (*o).funcpointer = Some(
                aubio_specdesc_decrease
                    as unsafe extern "C" fn(
                        _: *mut aubio_specdesc_t,
                        _: *const cvec_t,
                        _: *mut fvec_t,
                    ) -> (),
            )
        }
        15 => {
            (*o).funcpointer = Some(
                aubio_specdesc_rolloff
                    as unsafe extern "C" fn(
                        _: *mut aubio_specdesc_t,
                        _: *const cvec_t,
                        _: *mut fvec_t,
                    ) -> (),
            )
        }
        _ => {}
    }
    (*o).onset_type = onset_type;
    return o;
}
/* * deletion of a spectral descriptor

  \param o spectral descriptor object as returned by new_aubio_specdesc()

*/
#[no_mangle]
pub unsafe extern "C" fn del_aubio_specdesc(mut o: *mut aubio_specdesc_t) {
    match (*o).onset_type as u32 {
        3 => {
            del_fvec((*o).oldmag);
            del_fvec((*o).dev1);
            del_fvec((*o).theta1);
            del_fvec((*o).theta2);
        }
        4 | 5 => {
            del_fvec((*o).dev1);
            del_fvec((*o).theta1);
            del_fvec((*o).theta2);
            del_aubio_hist((*o).histog);
        }
        1 => {
            del_fvec((*o).oldmag);
            del_fvec((*o).dev1);
            del_aubio_hist((*o).histog);
        }
        6 | 7 | 8 => {
            del_fvec((*o).oldmag);
        }
        0 | 2 | _ => {}
    }
    free(o as *mut core::ffi::c_void);
}
