use crate::transpiled::specdesc::aubio_specdesc_t;

extern "C" {
    #[no_mangle]
    fn powf(_: f32, _: f32) -> f32;
    #[no_mangle]
    fn sqrtf(_: f32) -> f32;
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
#[derive(Copy, Clone)]
#[repr(C)]
pub struct cvec_t {
    pub length: uint_t,
    pub norm: *mut smpl_t,
    pub phas: *mut smpl_t,
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
#[no_mangle]
pub unsafe extern "C" fn cvec_sum(mut s: *const cvec_t) -> smpl_t {
    let mut j: uint_t = 0;
    let mut tmp: smpl_t = 0.0f64 as smpl_t;
    j = 0 as i32 as uint_t;
    while j < (*s).length {
        tmp += *(*s).norm.offset(j as isize);
        j = j.wrapping_add(1)
    }
    return tmp;
}
#[no_mangle]
pub unsafe extern "C" fn cvec_mean(mut s: *const cvec_t) -> smpl_t {
    return cvec_sum(s) / (*s).length as smpl_t;
}
#[no_mangle]
pub unsafe extern "C" fn cvec_centroid(mut spec: *const cvec_t) -> smpl_t {
    let mut sum: smpl_t = 0.0f64 as smpl_t;
    let mut sc: smpl_t = 0.0f64 as smpl_t;
    let mut j: uint_t = 0;
    sum = cvec_sum(spec);
    if sum as f64 == 0.0f64 {
        return 0.0f64 as smpl_t
    } else {
        j = 0 as i32 as uint_t;
        while j < (*spec).length {
            sc += j as smpl_t * *(*spec).norm.offset(j as isize);
            j = j.wrapping_add(1)
        }
        return sc / sum
    };
}
#[no_mangle]
pub unsafe extern "C" fn cvec_moment(mut spec: *const cvec_t,
                                     mut order: uint_t) -> smpl_t {
    let mut sum: smpl_t = 0.0f64 as smpl_t;
    let mut centroid: smpl_t = 0.0f64 as smpl_t;
    let mut sc: smpl_t = 0.0f64 as smpl_t;
    let mut j: uint_t = 0;
    sum = cvec_sum(spec);
    if sum as f64 == 0.0f64 {
        return 0.0f64 as smpl_t
    } else {
        centroid = cvec_centroid(spec);
        j = 0 as i32 as uint_t;
        while j < (*spec).length {
            sc +=
                powf(j as f32 - centroid, order as f32) *
                    *(*spec).norm.offset(j as isize);
            j = j.wrapping_add(1)
        }
        return sc / sum
    };
}
/*
  Copyright (C) 2007-2009 Paul Brossier <piem@aubio.org>

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
pub unsafe extern "C" fn aubio_specdesc_centroid(mut o: *mut aubio_specdesc_t,
                                                 mut spec: *const cvec_t,
                                                 mut desc: *mut fvec_t) {
    *(*desc).data.offset(0 as i32 as isize) = cvec_centroid(spec);
}
#[no_mangle]
pub unsafe extern "C" fn aubio_specdesc_spread(mut o: *mut aubio_specdesc_t,
                                               mut spec: *const cvec_t,
                                               mut desc: *mut fvec_t) {
    *(*desc).data.offset(0 as i32 as isize) =
        cvec_moment(spec, 2 as i32 as uint_t);
}
#[no_mangle]
pub unsafe extern "C" fn aubio_specdesc_skewness(mut o: *mut aubio_specdesc_t,
                                                 mut spec: *const cvec_t,
                                                 mut desc: *mut fvec_t) {
    let mut spread: smpl_t = 0.;
    spread = cvec_moment(spec, 2 as i32 as uint_t);
    if spread == 0 as i32 as f32 {
        *(*desc).data.offset(0 as i32 as isize) = 0.0f64 as smpl_t
    } else {
        *(*desc).data.offset(0 as i32 as isize) =
            cvec_moment(spec, 3 as i32 as uint_t);
        let ref mut fresh0 = *(*desc).data.offset(0 as i32 as isize);
        *fresh0 /= powf(sqrtf(spread), 3 as i32 as f32)
    };
}
#[no_mangle]
pub unsafe extern "C" fn aubio_specdesc_kurtosis(mut o: *mut aubio_specdesc_t,
                                                 mut spec: *const cvec_t,
                                                 mut desc: *mut fvec_t) {
    let mut spread: smpl_t = 0.;
    spread = cvec_moment(spec, 2 as i32 as uint_t);
    if spread == 0 as i32 as f32 {
        *(*desc).data.offset(0 as i32 as isize) = 0.0f64 as smpl_t
    } else {
        *(*desc).data.offset(0 as i32 as isize) =
            cvec_moment(spec, 4 as i32 as uint_t);
        let ref mut fresh1 = *(*desc).data.offset(0 as i32 as isize);
        *fresh1 /= spread * spread
    };
}
#[no_mangle]
pub unsafe extern "C" fn aubio_specdesc_slope(mut o: *mut aubio_specdesc_t,
                                              mut spec: *const cvec_t,
                                              mut desc: *mut fvec_t) {
    let mut j: uint_t = 0;
    let mut norm: smpl_t = 0 as i32 as smpl_t;
    let mut sum: smpl_t = 0.0f64 as smpl_t;
    // compute N * sum(j**2) - sum(j)**2
    j = 0 as i32 as uint_t;
    while j < (*spec).length {
        norm += j.wrapping_mul(j) as f32;
        j = j.wrapping_add(1)
    }
    norm *= (*spec).length as f32;
    // sum_0^N(j) = length * (length + 1) / 2
    norm =
        (norm as f64 -
             (*spec).length as f64 *
                 ((*spec).length as f64 - 1.0f64) / 2.0f64 *
                 ((*spec).length as f64 *
                      ((*spec).length as f64 - 1.0f64) / 2.0f64))
            as smpl_t;
    sum = cvec_sum(spec);
    *(*desc).data.offset(0 as i32 as isize) = 0.0f64 as smpl_t;
    if sum as f64 == 0.0f64 {
        return
    } else {
        j = 0 as i32 as uint_t;
        while j < (*spec).length {
            let ref mut fresh2 =
                *(*desc).data.offset(0 as i32 as isize);
            *fresh2 += j as f32 * *(*spec).norm.offset(j as isize);
            j = j.wrapping_add(1)
        }
        let ref mut fresh3 = *(*desc).data.offset(0 as i32 as isize);
        *fresh3 *= (*spec).length as f32;
        let ref mut fresh4 = *(*desc).data.offset(0 as i32 as isize);
        *fresh4 =
            (*fresh4 as f64 -
                 (sum * (*spec).length as f32 *
                      (*spec).length.wrapping_sub(1 as i32 as
                                                      u32) as
                          f32) as f64 / 2.0f64) as
                smpl_t;
        let ref mut fresh5 = *(*desc).data.offset(0 as i32 as isize);
        *fresh5 /= norm;
        let ref mut fresh6 = *(*desc).data.offset(0 as i32 as isize);
        *fresh6 /= sum
    };
}
#[no_mangle]
pub unsafe extern "C" fn aubio_specdesc_decrease(mut o: *mut aubio_specdesc_t,
                                                 mut spec: *const cvec_t,
                                                 mut desc: *mut fvec_t) {
    let mut j: uint_t = 0;
    let mut sum: smpl_t = 0.;
    sum = cvec_sum(spec);
    *(*desc).data.offset(0 as i32 as isize) =
        0 as i32 as smpl_t;
    if sum as f64 == 0.0f64 {
        return
    } else {
        sum -= *(*spec).norm.offset(0 as i32 as isize);
        j = 1 as i32 as uint_t;
        while j < (*spec).length {
            let ref mut fresh7 =
                *(*desc).data.offset(0 as i32 as isize);
            *fresh7 +=
                (*(*spec).norm.offset(j as isize) -
                     *(*spec).norm.offset(0 as i32 as isize)) /
                    j as f32;
            j = j.wrapping_add(1)
        }
        let ref mut fresh8 = *(*desc).data.offset(0 as i32 as isize);
        *fresh8 /= sum
    };
}
#[no_mangle]
pub unsafe extern "C" fn aubio_specdesc_rolloff(mut o: *mut aubio_specdesc_t,
                                                mut spec: *const cvec_t,
                                                mut desc: *mut fvec_t) {
    let mut j: uint_t = 0;
    let mut cumsum: smpl_t = 0.;
    let mut rollsum: smpl_t = 0.;
    cumsum = 0.0f64 as smpl_t;
    rollsum = 0.0f64 as smpl_t;
    j = 0 as i32 as uint_t;
    while j < (*spec).length {
        cumsum +=
            *(*spec).norm.offset(j as isize) *
                *(*spec).norm.offset(j as isize);
        j = j.wrapping_add(1)
    }
    if cumsum == 0 as i32 as f32 {
        *(*desc).data.offset(0 as i32 as isize) = 0.0f64 as smpl_t
    } else {
        cumsum = (cumsum as f64 * 0.95f64) as smpl_t;
        j = 0 as i32 as uint_t;
        while rollsum < cumsum {
            rollsum +=
                *(*spec).norm.offset(j as isize) *
                    *(*spec).norm.offset(j as isize);
            j = j.wrapping_add(1)
        }
        *(*desc).data.offset(0 as i32 as isize) = j as smpl_t
    };
}
