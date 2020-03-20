use ::libc;
extern "C" {
    pub type _aubio_pvoc_t;
    pub type _aubio_filter_t;
    pub type _aubio_pitchmcomb_t;
    pub type _aubio_pitchyin_t;
    pub type _aubio_pitchfcomb_t;
    pub type _aubio_pitchschmitt_t;
    pub type _aubio_pitchyinfft_t;
    pub type _aubio_pitchyinfast_t;
    pub type _aubio_pitchspecacf_t;
    #[no_mangle]
    fn free(_: *mut libc::c_void);
    #[no_mangle]
    fn calloc(_: libc::c_ulong, _: libc::c_ulong) -> *mut libc::c_void;
    #[no_mangle]
    fn strcmp(_: *const libc::c_char, _: *const libc::c_char) -> libc::c_int;
    /* file interface */
    /* strings */
    /* Error reporting */
    /* Logging */
    /* * internal logging function, defined in utils/log.c */
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
    /* * convert frequency bin to frequency (Hz) */
    #[no_mangle]
    fn aubio_bintofreq(bin: smpl_t, samplerate: smpl_t, fftsize: smpl_t)
     -> smpl_t;
    /* * convert frequency (Hz) to frequency bin */
    #[no_mangle]
    fn aubio_freqtobin(freq: smpl_t, samplerate: smpl_t, fftsize: smpl_t)
     -> smpl_t;
    /* * convert frequency (Hz) to midi value (0-128) */
    #[no_mangle]
    fn aubio_freqtomidi(freq: smpl_t) -> smpl_t;
    /* * check if buffer level in dB SPL is under a given threshold

  \param v vector to get level from
  \param threshold threshold in dB SPL

  \return 0 if level is under the given threshold, 1 otherwise

*/
    #[no_mangle]
    fn aubio_silence_detection(v: *const fvec_t, threshold: smpl_t) -> uint_t;
    /* * cvec_t buffer creation function

  This function creates a cvec_t structure holding two arrays of size
  [length/2+1], corresponding to the norm and phase values of the
  spectral frame. The length stored in the structure is the actual size of both
  arrays, not the length of the complex and symmetrical vector, specified as
  creation argument.

  \param length the length of the buffer to create

*/
    #[no_mangle]
    fn new_cvec(length: uint_t) -> *mut cvec_t;
    /* * cvec_t buffer deletion function

  \param s buffer to delete as returned by new_cvec()

*/
    #[no_mangle]
    fn del_cvec(s: *mut cvec_t);
    /* * create phase vocoder object

  \param win_s size of analysis buffer (and length the FFT transform)
  \param hop_s step size between two consecutive analysis

*/
    #[no_mangle]
    fn new_aubio_pvoc(win_s: uint_t, hop_s: uint_t) -> *mut aubio_pvoc_t;
    /* * delete phase vocoder object

  \param pv phase vocoder object as returned by new_aubio_pvoc

*/
    #[no_mangle]
    fn del_aubio_pvoc(pv: *mut aubio_pvoc_t);
    /* * compute spectral frame

  This function accepts an input vector of size [hop_s]. The
  analysis buffer is rotated and filled with the new data. After windowing of
  this signal window, the Fourier transform is computed and returned in
  fftgrain as two vectors, magnitude and phase.

  \param pv phase vocoder object as returned by new_aubio_pvoc
  \param in new input signal (hop_s long)
  \param fftgrain output spectral frame

*/
    #[no_mangle]
    fn aubio_pvoc_do(pv: *mut aubio_pvoc_t, in_0: *const fvec_t,
                     fftgrain: *mut cvec_t);
    /* * filter input vector (out-of-place)

  \param f filter object as returned by new_aubio_filter()
  \param in input vector to filter
  \param out output vector to store filtered input

*/
    #[no_mangle]
    fn aubio_filter_do_outplace(f: *mut aubio_filter_t, in_0: *const fvec_t,
                                out: *mut fvec_t);
    /* * delete a filter object

  \param f filter object to delete

*/
    #[no_mangle]
    fn del_aubio_filter(f: *mut aubio_filter_t);
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

  C-weighting filter coefficients

  This file creates a C-weighting digital filter, which reduces low and high
  frequencies and enhance the middle ones to reflect the ability of the human
  hearing.

  The implementation is based on the following standard:

    - IEC/CD 1672: Electroacoustics-Sound Level Meters, IEC, Geneva, Nov.  1996,
  for A- and C-weighting filters.

  See also:

    - <a href="http://en.wikipedia.org/wiki/A-weighting">A-Weighting on
  Wikipedia</a>
    - <a href="http://en.wikipedia.org/wiki/Weighting_filter">Weighting filter on
  Wikipedia</a>
    - <a href="http://www.mathworks.com/matlabcentral/fileexchange/69">Christophe
  Couvreur's 'octave' toolbox</a>

  The coefficients in this file have been computed using Christophe Couvreur's
  scripts in octave 3.0 (debian package 1:3.0.5-6+b2 with octave-signal
  1.0.9-1+b1 on i386), with <pre> [b, a] = cdsign(1/Fs) </pre> for various
  sampling frequencies (8000, 11025, 16000, 22050, 24000, 32000, 44100, 48000,
  88200, 96000, and 192000 Hz).

  The sampling frequency should normally be higher than 20kHz, but most common
  file sampling rates have been included for completeness.

  \example temporal/test-c_weighting.c

*/
    /* * create new C-design filter

  \param samplerate sampling frequency of the signal to filter. Should be one of
  8000, 11025, 16000, 22050, 24000, 32000, 44100, 48000, 88200, 96000, and
  192000 Hz

  \return a new filter object

*/
    #[no_mangle]
    fn new_aubio_filter_c_weighting(samplerate: uint_t)
     -> *mut aubio_filter_t;
    /* * execute pitch detection on an input spectral frame

  \param p pitch detection object as returned by new_aubio_pitchmcomb
  \param in_fftgrain input signal spectrum as computed by aubio_pvoc_do
  \param out_cands pitch candidate frequenciess, in bins

*/
    #[no_mangle]
    fn aubio_pitchmcomb_do(p: *mut aubio_pitchmcomb_t,
                           in_fftgrain: *const cvec_t,
                           out_cands: *mut fvec_t);
    /* * creation of the pitch detection object

  \param buf_size size of the input buffer to analyse
  \param hop_size step size between two consecutive analysis instant

*/
    #[no_mangle]
    fn new_aubio_pitchmcomb(buf_size: uint_t, hop_size: uint_t)
     -> *mut aubio_pitchmcomb_t;
    /* * deletion of the pitch detection object

  \param p pitch detection object as returned by new_aubio_pitchfcomb

*/
    #[no_mangle]
    fn del_aubio_pitchmcomb(p: *mut aubio_pitchmcomb_t);
    /* * creation of the pitch detection object

  \param buf_size size of the input buffer to analyse

*/
    #[no_mangle]
    fn new_aubio_pitchyin(buf_size: uint_t) -> *mut aubio_pitchyin_t;
    /* * deletion of the pitch detection object

  \param o pitch detection object as returned by new_aubio_pitchyin()

*/
    #[no_mangle]
    fn del_aubio_pitchyin(o: *mut aubio_pitchyin_t);
    /* * execute pitch detection an input buffer

  \param o pitch detection object as returned by new_aubio_pitchyin()
  \param samples_in input signal vector (length as specified at creation time)
  \param cands_out pitch period candidates, in samples

*/
    #[no_mangle]
    fn aubio_pitchyin_do(o: *mut aubio_pitchyin_t, samples_in: *const fvec_t,
                         cands_out: *mut fvec_t);
    /* * set tolerance parameter for YIN algorithm

  \param o YIN pitch detection object
  \param tol tolerance parameter for minima selection [default 0.15]

*/
    #[no_mangle]
    fn aubio_pitchyin_set_tolerance(o: *mut aubio_pitchyin_t, tol: smpl_t)
     -> uint_t;
    /* * get tolerance parameter for YIN algorithm

  \param o YIN pitch detection object
  \return tolerance parameter for minima selection [default 0.15]

*/
    #[no_mangle]
    fn aubio_pitchyin_get_tolerance(o: *mut aubio_pitchyin_t) -> smpl_t;
    /* * get current confidence of YIN algorithm

  \param o YIN pitch detection object
  \return confidence parameter

*/
    #[no_mangle]
    fn aubio_pitchyin_get_confidence(o: *mut aubio_pitchyin_t) -> smpl_t;
    /* * execute pitch detection on an input buffer

  \param p pitch detection object as returned by new_aubio_pitchfcomb
  \param input input signal window (length as specified at creation time)
  \param output pitch candidates in bins

*/
    #[no_mangle]
    fn aubio_pitchfcomb_do(p: *mut aubio_pitchfcomb_t, input: *const fvec_t,
                           output: *mut fvec_t);
    /* * creation of the pitch detection object

  \param buf_size size of the input buffer to analyse
  \param hop_size step size between two consecutive analysis instant

*/
    #[no_mangle]
    fn new_aubio_pitchfcomb(buf_size: uint_t, hop_size: uint_t)
     -> *mut aubio_pitchfcomb_t;
    /* * deletion of the pitch detection object

  \param p pitch detection object as returned by new_aubio_pitchfcomb

*/
    #[no_mangle]
    fn del_aubio_pitchfcomb(p: *mut aubio_pitchfcomb_t);
    /* * execute pitch detection on an input buffer

  \param p pitch detection object as returned by new_aubio_pitchschmitt
  \param samples_in input signal vector (length as specified at creation time)
  \param cands_out pitch period estimates, in samples

*/
    #[no_mangle]
    fn aubio_pitchschmitt_do(p: *mut aubio_pitchschmitt_t,
                             samples_in: *const fvec_t,
                             cands_out: *mut fvec_t);
    /* * creation of the pitch detection object

  \param buf_size size of the input buffer to analyse

*/
    #[no_mangle]
    fn new_aubio_pitchschmitt(buf_size: uint_t) -> *mut aubio_pitchschmitt_t;
    /* * deletion of the pitch detection object

  \param p pitch detection object as returned by new_aubio_pitchschmitt

*/
    #[no_mangle]
    fn del_aubio_pitchschmitt(p: *mut aubio_pitchschmitt_t);
    /* * execute pitch detection on an input buffer

  \param o pitch detection object as returned by new_aubio_pitchyinfft
  \param samples_in input signal vector (length as specified at creation time)
  \param cands_out pitch period candidates, in samples

*/
    #[no_mangle]
    fn aubio_pitchyinfft_do(o: *mut aubio_pitchyinfft_t,
                            samples_in: *const fvec_t,
                            cands_out: *mut fvec_t);
    /* * creation of the pitch detection object

  \param samplerate samplerate of the input signal
  \param buf_size size of the input buffer to analyse

*/
    #[no_mangle]
    fn new_aubio_pitchyinfft(samplerate: uint_t, buf_size: uint_t)
     -> *mut aubio_pitchyinfft_t;
    /* * deletion of the pitch detection object

  \param o pitch detection object as returned by new_aubio_pitchyinfft()

*/
    #[no_mangle]
    fn del_aubio_pitchyinfft(o: *mut aubio_pitchyinfft_t);
    /* * get tolerance parameter for YIN algorithm

  \param o YIN pitch detection object

  \return tolerance parameter for minima selection [default 0.15]

*/
    #[no_mangle]
    fn aubio_pitchyinfft_get_tolerance(o: *mut aubio_pitchyinfft_t) -> smpl_t;
    /* * set tolerance parameter for YIN algorithm

  \param o YIN pitch detection object
  \param tol tolerance parameter for minima selection [default 0.15]

*/
    #[no_mangle]
    fn aubio_pitchyinfft_set_tolerance(o: *mut aubio_pitchyinfft_t,
                                       tol: smpl_t) -> uint_t;
    /* * get current confidence of YIN algorithm

  \param o YIN pitch detection object
  \return confidence parameter

*/
    #[no_mangle]
    fn aubio_pitchyinfft_get_confidence(o: *mut aubio_pitchyinfft_t)
     -> smpl_t;
    /* * creation of the pitch detection object

  \param buf_size size of the input buffer to analyse

*/
    #[no_mangle]
    fn new_aubio_pitchyinfast(buf_size: uint_t) -> *mut aubio_pitchyinfast_t;
    /* * deletion of the pitch detection object

  \param o pitch detection object as returned by new_aubio_pitchyin()

*/
    #[no_mangle]
    fn del_aubio_pitchyinfast(o: *mut aubio_pitchyinfast_t);
    /* * execute pitch detection an input buffer

  \param o pitch detection object as returned by new_aubio_pitchyin()
  \param samples_in input signal vector (length as specified at creation time)
  \param cands_out pitch period candidates, in samples

*/
    #[no_mangle]
    fn aubio_pitchyinfast_do(o: *mut aubio_pitchyinfast_t,
                             samples_in: *const fvec_t,
                             cands_out: *mut fvec_t);
    /* * set tolerance parameter for YIN algorithm

  \param o YIN pitch detection object
  \param tol tolerance parameter for minima selection [default 0.15]

*/
    #[no_mangle]
    fn aubio_pitchyinfast_set_tolerance(o: *mut aubio_pitchyinfast_t,
                                        tol: smpl_t) -> uint_t;
    /* * get tolerance parameter for YIN algorithm

  \param o YIN pitch detection object
  \return tolerance parameter for minima selection [default 0.15]

*/
    #[no_mangle]
    fn aubio_pitchyinfast_get_tolerance(o: *mut aubio_pitchyinfast_t)
     -> smpl_t;
    /* * get current confidence of YIN algorithm

  \param o YIN pitch detection object
  \return confidence parameter

*/
    #[no_mangle]
    fn aubio_pitchyinfast_get_confidence(o: *mut aubio_pitchyinfast_t)
     -> smpl_t;
    /* * execute pitch detection on an input buffer

  \param o pitch detection object as returned by new_aubio_pitchspecacf
  \param samples_in input signal vector (length as specified at creation time)
  \param cands_out pitch period candidates, in samples

*/
    #[no_mangle]
    fn aubio_pitchspecacf_do(o: *mut aubio_pitchspecacf_t,
                             samples_in: *const fvec_t,
                             cands_out: *mut fvec_t);
    /* * creation of the pitch detection object

  \param buf_size size of the input buffer to analyse

*/
    #[no_mangle]
    fn new_aubio_pitchspecacf(buf_size: uint_t) -> *mut aubio_pitchspecacf_t;
    /* * deletion of the pitch detection object

  \param o pitch detection object as returned by new_aubio_pitchspecacf()

*/
    #[no_mangle]
    fn del_aubio_pitchspecacf(o: *mut aubio_pitchspecacf_t);
    /* * get tolerance parameter for `specacf` pitch detection object

  \param o pitch detection object

  \return tolerance parameter for minima selection [default 1.]

*/
    #[no_mangle]
    fn aubio_pitchspecacf_get_tolerance(o: *const aubio_pitchspecacf_t)
     -> smpl_t;
    /* * set tolerance parameter for `specacf` pitch detection object

  \param o pitch detection object
  \param tol tolerance parameter for minima selection [default 1.]

  \return `1` on error, `0` on success

*/
    #[no_mangle]
    fn aubio_pitchspecacf_set_tolerance(o: *mut aubio_pitchspecacf_t,
                                        tol: smpl_t) -> uint_t;
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

  Phase vocoder object

  This object implements a phase vocoder. The spectral frames are computed
  using a HanningZ window and a swapped version of the signal to simplify the
  phase relationships across frames. The window sizes and overlap are specified
  at creation time.

  \example spectral/test-phasevoc.c

*/
/* * phasevocoder object */
pub type aubio_pvoc_t = _aubio_pvoc_t;
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

  Digital filter

  This object stores a digital filter of order \f$n\f$.
  It contains the following data:
    - \f$ n*1 b_i \f$ feedforward coefficients
    - \f$ n*1 a_i \f$ feedback coefficients
    - \f$ n*c x_i \f$ input signal
    - \f$ n*c y_i \f$ output signal

  For convenience, the samplerate of the input signal is also stored in the
  object.

  Feedforward and feedback parameters can be modified using
  aubio_filter_get_feedback() and aubio_filter_get_feedforward().

  The function aubio_filter_do_outplace() computes the following output signal
  \f$ y[n] \f$ from the input signal \f$ x[n] \f$:

  \f{eqnarray*}{
     y[n] = b_0 x[n] & + & b_1 x[n-1] + b_2 x[n-2] + ... + b_P x[n-P] \\
                     & - & a_1 y[n-1] - a_2 y[n-2] - ... - a_P y[n-P] \\
  \f}

  The function aubio_filter_do() executes the same computation but modifies
  directly the input signal (in-place).

  The function aubio_filter_do_filtfilt() version runs the filter twice, first
  forward then backward, to compensate with the phase shifting of the forward
  operation.

  Some convenience functions are provided:
    - new_aubio_filter_a_weighting() and aubio_filter_set_a_weighting(),
    - new_aubio_filter_c_weighting() and aubio_filter_set_c_weighting().
    - new_aubio_filter_biquad() and aubio_filter_set_biquad().

  \example temporal/test-filter.c

*/
/* * Digital filter

*/
pub type aubio_filter_t = _aubio_filter_t;
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

  Pitch detection using multiple-comb filter

  This fundamental frequency estimation algorithm implements spectral
  flattening, multi-comb filtering and peak histogramming.

  This method was designed by Juan P. Bello and described in:

  Juan-Pablo Bello. ``Towards the Automated Analysis of Simple Polyphonic
  Music''.  PhD thesis, Centre for Digital Music, Queen Mary University of
  London, London, UK, 2003.

  \example pitch/test-pitchmcomb.c

*/
/* * pitch detection object */
pub type aubio_pitchmcomb_t = _aubio_pitchmcomb_t;
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

  Pitch detection using the YIN algorithm

  This algorithm was developed by A. de Cheveigne and H. Kawahara and
  published in:

  De Cheveigné, A., Kawahara, H. (2002) "YIN, a fundamental frequency
  estimator for speech and music", J. Acoust. Soc. Am. 111, 1917-1930.

  see http://recherche.ircam.fr/equipes/pcm/pub/people/cheveign.html
      http://recherche.ircam.fr/equipes/pcm/cheveign/ps/2002_JASA_YIN_proof.pdf

  \example pitch/test-pitchyin.c

*/
/* * pitch detection object */
pub type aubio_pitchyin_t = _aubio_pitchyin_t;
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

   Pitch detection using a fast harmonic comb filter

   This pitch extraction method implements a fast harmonic comb filter to
   determine the fundamental frequency of a harmonic sound.

   This file was derived from the tuneit project, written by Mario Lang to
   detect the fundamental frequency of a sound.

   See http://delysid.org/tuneit.html

   \example pitch/test-pitchfcomb.c

*/
/* * pitch detection object */
pub type aubio_pitchfcomb_t = _aubio_pitchfcomb_t;
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

   Pitch detection using a Schmitt trigger

   This pitch extraction method implements a Schmitt trigger to estimate the
   period of a signal.

   This file was derived from the tuneit project, written by Mario Lang to
   detect the fundamental frequency of a sound.

   See http://delysid.org/tuneit.html

   \example pitch/test-pitchschmitt.c

*/
/* * pitch detection object */
pub type aubio_pitchschmitt_t = _aubio_pitchschmitt_t;
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

  Pitch detection using a spectral implementation of the YIN algorithm

  This algorithm was derived from the YIN algorithm. In this implementation, a
  Fourier transform is used to compute a tapered square difference function,
  which allows spectral weighting. Because the difference function is tapered,
  the selection of the period is simplified.

  Paul Brossier, [Automatic annotation of musical audio for interactive
  systems](http://aubio.org/phd/), Chapter 3, Pitch Analysis, PhD thesis,
  Centre for Digital music, Queen Mary University of London, London, UK, 2006.

  \example pitch/test-pitchyinfft.c

*/
/* * pitch detection object */
pub type aubio_pitchyinfft_t = _aubio_pitchyinfft_t;
/*
  Copyright (C) 2003-2017 Paul Brossier <piem@aubio.org>

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

  Pitch detection using YIN algorithm (fast implementation)

  This algorithm was developed by A. de Cheveigne and H. Kawahara and
  published in:

  De Cheveigné, A., Kawahara, H. (2002) "YIN, a fundamental frequency
  estimator for speech and music", J. Acoust. Soc. Am. 111, 1917-1930.

  This implementation compute the autocorrelation function using time domain
  convolution computed in the spectral domain.

  see http://recherche.ircam.fr/equipes/pcm/pub/people/cheveign.html
      http://recherche.ircam.fr/equipes/pcm/cheveign/ps/2002_JASA_YIN_proof.pdf

*/
/* * pitch detection object */
pub type aubio_pitchyinfast_t = _aubio_pitchyinfast_t;
/*
  Copyright (C) 2013 Paul Brossier <piem@aubio.org>

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

  Pitch detection using spectral auto correlation

  This algorithm implements pitch detection by computing the autocorrelation
  function as the cosine transform of the square spectral magnitudes.

  Anssi Klapuri. Qualitative and quantitative aspects in the design of
  periodicity esti- mation algorithms. In Proceedings of the European Signal
  Processing Conference (EUSIPCO), 2000.

  Paul Brossier, [Automatic annotation of musical audio for interactive
  systems](http://aubio.org/phd/), Chapter 3, Pitch Analysis, Autocorrelation,
  pp. 75-77, PhD thesis, Centre for Digital music, Queen Mary University of
  London, London, UK, 2006.

  \example pitch/test-pitchspecacf.c

*/
/* * pitch detection object */
pub type aubio_pitchspecacf_t = _aubio_pitchspecacf_t;
/* * generic pitch detection structure */
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _aubio_pitch_t {
    pub type_0: aubio_pitch_type,
    pub mode: aubio_pitch_mode,
    pub samplerate: uint_t,
    pub bufsize: uint_t,
    pub p_object: *mut libc::c_void,
    pub filter: *mut aubio_filter_t,
    pub filtered: *mut fvec_t,
    pub pv: *mut aubio_pvoc_t,
    pub fftgrain: *mut cvec_t,
    pub buf: *mut fvec_t,
    pub detect_cb: aubio_pitch_detect_t,
    pub conv_cb: aubio_pitch_convert_t,
    pub conf_cb: aubio_pitch_get_conf_t,
    pub silence: smpl_t,
}
/* * callback to fetch the confidence of the algorithm */
pub type aubio_pitch_get_conf_t
    =
    Option<unsafe extern "C" fn(_: *mut libc::c_void) -> smpl_t>;
/* * callback to convert pitch from one unit to another, defined below */
pub type aubio_pitch_convert_t
    =
    Option<unsafe extern "C" fn(_: smpl_t, _: uint_t, _: uint_t) -> smpl_t>;
/* * callback to get pitch candidate, defined below */
pub type aubio_pitch_detect_t
    =
    Option<unsafe extern "C" fn(_: *mut aubio_pitch_t, _: *const fvec_t,
                                _: *mut fvec_t) -> ()>;
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

  Pitch detection object

  This file creates the objects required for the computation of the selected
  pitch detection algorithm and output the results, in midi note or Hz.

  \section pitch Pitch detection methods

  A list of the pitch detection methods currently available follows.

  \b \p default : use the default method

  Currently, the default method is set to \p yinfft .

  \b \p schmitt : Schmitt trigger

  This pitch extraction method implements a Schmitt trigger to estimate the
  period of a signal.

  This file was derived from the tuneit project, written by Mario Lang to
  detect the fundamental frequency of a sound.

  See http://delysid.org/tuneit.html

  \b \p fcomb : a fast harmonic comb filter

  This pitch extraction method implements a fast harmonic comb filter to
  determine the fundamental frequency of a harmonic sound.

  This file was derived from the tuneit project, written by Mario Lang to
  detect the fundamental frequency of a sound.

  See http://delysid.org/tuneit.html

  \b \p mcomb : multiple-comb filter

  This fundamental frequency estimation algorithm implements spectral
  flattening, multi-comb filtering and peak histogramming.

  This method was designed by Juan P. Bello and described in:

  Juan-Pablo Bello. ``Towards the Automated Analysis of Simple Polyphonic
  Music''.  PhD thesis, Centre for Digital Music, Queen Mary University of
  London, London, UK, 2003.

  \b \p yin : YIN algorithm

  This algorithm was developed by A. de Cheveigne and H. Kawahara and
  published in:

  De Cheveigné, A., Kawahara, H. (2002) "YIN, a fundamental frequency
  estimator for speech and music", J. Acoust. Soc. Am. 111, 1917-1930.

  see http://recherche.ircam.fr/equipes/pcm/pub/people/cheveign.html

  \b \p yinfast : Yinfast algorithm

  This algorithm is equivalent to the YIN algorithm, but computed in the
  spectral domain for efficiency. See also `python/demos/demo_yin_compare.py`.

  \b \p yinfft : Yinfft algorithm

  This algorithm was derived from the YIN algorithm. In this implementation, a
  Fourier transform is used to compute a tapered square difference function,
  which allows spectral weighting. Because the difference function is tapered,
  the selection of the period is simplified.

  Paul Brossier, [Automatic annotation of musical audio for interactive
  systems](http://aubio.org/phd/), Chapter 3, Pitch Analysis, PhD thesis,
  Centre for Digital music, Queen Mary University of London, London, UK, 2006.

  \example pitch/test-pitch.c
  \example examples/aubiopitch.c

*/
/* * pitch detection object */
pub type aubio_pitch_t = _aubio_pitch_t;
pub type aubio_pitch_mode = libc::c_uint;
pub const aubio_pitchm_default: aubio_pitch_mode = 0;
pub const aubio_pitchm_bin: aubio_pitch_mode = 3;
pub const aubio_pitchm_cent: aubio_pitch_mode = 2;
pub const aubio_pitchm_midi: aubio_pitch_mode = 1;
pub const aubio_pitchm_freq: aubio_pitch_mode = 0;
pub type aubio_pitch_type = libc::c_uint;
pub const aubio_pitcht_default: aubio_pitch_type = 4;
pub const aubio_pitcht_specacf: aubio_pitch_type = 6;
pub const aubio_pitcht_yinfast: aubio_pitch_type = 5;
pub const aubio_pitcht_yinfft: aubio_pitch_type = 4;
pub const aubio_pitcht_fcomb: aubio_pitch_type = 3;
pub const aubio_pitcht_schmitt: aubio_pitch_type = 2;
pub const aubio_pitcht_mcomb: aubio_pitch_type = 1;
pub const aubio_pitcht_yin: aubio_pitch_type = 0;
/* * creation of the pitch detection object

  \param method set pitch detection algorithm
  \param buf_size size of the input buffer to analyse
  \param hop_size step size between two consecutive analysis instant
  \param samplerate sampling rate of the signal

  \return newly created ::aubio_pitch_t

*/
#[no_mangle]
pub unsafe extern "C" fn new_aubio_pitch(mut pitch_mode: *const char_t,
                                         mut bufsize: uint_t,
                                         mut hopsize: uint_t,
                                         mut samplerate: uint_t)
 -> *mut aubio_pitch_t {
    let mut current_block: u64;
    let mut p: *mut aubio_pitch_t =
        calloc(::std::mem::size_of::<aubio_pitch_t>() as libc::c_ulong,
               1 as libc::c_int as libc::c_ulong) as *mut aubio_pitch_t;
    let mut pitch_type: aubio_pitch_type = aubio_pitcht_yin;
    if pitch_mode.is_null() {
        aubio_log(AUBIO_LOG_ERR as libc::c_int,
                  b"AUBIO ERROR: pitch: can not use \xe2\x80\x98NULL\xe2\x80\x98 for pitch detection method\n\x00"
                      as *const u8 as *const libc::c_char);
    } else {
        if strcmp(pitch_mode,
                  b"mcomb\x00" as *const u8 as *const libc::c_char) ==
               0 as libc::c_int {
            pitch_type = aubio_pitcht_mcomb;
            current_block = 5689001924483802034;
        } else if strcmp(pitch_mode,
                         b"yinfast\x00" as *const u8 as *const libc::c_char)
                      == 0 as libc::c_int {
            pitch_type = aubio_pitcht_yinfast;
            current_block = 5689001924483802034;
        } else if strcmp(pitch_mode,
                         b"yinfft\x00" as *const u8 as *const libc::c_char) ==
                      0 as libc::c_int {
            pitch_type = aubio_pitcht_yinfft;
            current_block = 5689001924483802034;
        } else if strcmp(pitch_mode,
                         b"yin\x00" as *const u8 as *const libc::c_char) ==
                      0 as libc::c_int {
            pitch_type = aubio_pitcht_yin;
            current_block = 5689001924483802034;
        } else if strcmp(pitch_mode,
                         b"schmitt\x00" as *const u8 as *const libc::c_char)
                      == 0 as libc::c_int {
            pitch_type = aubio_pitcht_schmitt;
            current_block = 5689001924483802034;
        } else if strcmp(pitch_mode,
                         b"fcomb\x00" as *const u8 as *const libc::c_char) ==
                      0 as libc::c_int {
            pitch_type = aubio_pitcht_fcomb;
            current_block = 5689001924483802034;
        } else if strcmp(pitch_mode,
                         b"specacf\x00" as *const u8 as *const libc::c_char)
                      == 0 as libc::c_int {
            pitch_type = aubio_pitcht_specacf;
            current_block = 5689001924483802034;
        } else if strcmp(pitch_mode,
                         b"default\x00" as *const u8 as *const libc::c_char)
                      == 0 as libc::c_int {
            pitch_type = aubio_pitcht_default;
            current_block = 5689001924483802034;
        } else {
            aubio_log(AUBIO_LOG_ERR as libc::c_int,
                      b"AUBIO ERROR: pitch: unknown pitch detection method \xe2\x80\x98%s\xe2\x80\x99\n\x00"
                          as *const u8 as *const libc::c_char, pitch_mode);
            current_block = 9082146709813991249;
        }
        match current_block {
            9082146709813991249 => { }
            _ =>
            // check parameters are valid
            {
                if (hopsize as sint_t) < 1 as libc::c_int {
                    aubio_log(AUBIO_LOG_ERR as libc::c_int,
                              b"AUBIO ERROR: pitch: got hopsize %d, but can not be < 1\n\x00"
                                  as *const u8 as *const libc::c_char,
                              hopsize);
                } else if (bufsize as sint_t) < 1 as libc::c_int {
                    aubio_log(AUBIO_LOG_ERR as libc::c_int,
                              b"AUBIO ERROR: pitch: got buffer_size %d, but can not be < 1\n\x00"
                                  as *const u8 as *const libc::c_char,
                              bufsize);
                } else if bufsize < hopsize {
                    aubio_log(AUBIO_LOG_ERR as libc::c_int,
                              b"AUBIO ERROR: pitch: hop size (%d) is larger than win size (%d)\n\x00"
                                  as *const u8 as *const libc::c_char,
                              hopsize, bufsize);
                } else if (samplerate as sint_t) < 1 as libc::c_int {
                    aubio_log(AUBIO_LOG_ERR as libc::c_int,
                              b"AUBIO ERROR: pitch: samplerate (%d) can not be < 1\n\x00"
                                  as *const u8 as *const libc::c_char,
                              samplerate);
                } else {
                    (*p).samplerate = samplerate;
                    (*p).type_0 = pitch_type;
                    aubio_pitch_set_unit(p,
                                         b"default\x00" as *const u8 as
                                             *const libc::c_char);
                    (*p).bufsize = bufsize;
                    (*p).silence = -50.0f64 as smpl_t;
                    (*p).conf_cb = None;
                    match (*p).type_0 as libc::c_uint {
                        0 => {
                            (*p).buf = new_fvec(bufsize);
                            (*p).p_object =
                                new_aubio_pitchyin(bufsize) as
                                    *mut libc::c_void;
                            if (*p).p_object.is_null() {
                                current_block = 9082146709813991249;
                            } else {
                                (*p).detect_cb =
                                    Some(aubio_pitch_do_yin as
                                             unsafe extern "C" fn(_:
                                                                      *mut aubio_pitch_t,
                                                                  _:
                                                                      *const fvec_t,
                                                                  _:
                                                                      *mut fvec_t)
                                                 -> ());
                                (*p).conf_cb =
                                    ::std::mem::transmute::<Option<unsafe extern "C" fn(_:
                                                                                            *mut aubio_pitchyin_t)
                                                                       ->
                                                                           smpl_t>,
                                                            aubio_pitch_get_conf_t>(Some(aubio_pitchyin_get_confidence
                                                                                             as
                                                                                             unsafe extern "C" fn(_:
                                                                                                                      *mut aubio_pitchyin_t)
                                                                                                 ->
                                                                                                     smpl_t));
                                aubio_pitchyin_set_tolerance((*p).p_object as
                                                                 *mut aubio_pitchyin_t,
                                                             0.15f64 as
                                                                 smpl_t);
                                current_block = 9437375157805982253;
                            }
                        }
                        1 => {
                            (*p).filtered = new_fvec(hopsize);
                            (*p).pv = new_aubio_pvoc(bufsize, hopsize);
                            if (*p).pv.is_null() {
                                current_block = 9082146709813991249;
                            } else {
                                (*p).fftgrain = new_cvec(bufsize);
                                (*p).p_object =
                                    new_aubio_pitchmcomb(bufsize, hopsize) as
                                        *mut libc::c_void;
                                (*p).filter =
                                    new_aubio_filter_c_weighting(samplerate);
                                (*p).detect_cb =
                                    Some(aubio_pitch_do_mcomb as
                                             unsafe extern "C" fn(_:
                                                                      *mut aubio_pitch_t,
                                                                  _:
                                                                      *const fvec_t,
                                                                  _:
                                                                      *mut fvec_t)
                                                 -> ());
                                current_block = 9437375157805982253;
                            }
                        }
                        3 => {
                            (*p).buf = new_fvec(bufsize);
                            (*p).p_object =
                                new_aubio_pitchfcomb(bufsize, hopsize) as
                                    *mut libc::c_void;
                            if (*p).p_object.is_null() {
                                current_block = 9082146709813991249;
                            } else {
                                (*p).detect_cb =
                                    Some(aubio_pitch_do_fcomb as
                                             unsafe extern "C" fn(_:
                                                                      *mut aubio_pitch_t,
                                                                  _:
                                                                      *const fvec_t,
                                                                  _:
                                                                      *mut fvec_t)
                                                 -> ());
                                current_block = 9437375157805982253;
                            }
                        }
                        2 => {
                            (*p).buf = new_fvec(bufsize);
                            (*p).p_object =
                                new_aubio_pitchschmitt(bufsize) as
                                    *mut libc::c_void;
                            (*p).detect_cb =
                                Some(aubio_pitch_do_schmitt as
                                         unsafe extern "C" fn(_:
                                                                  *mut aubio_pitch_t,
                                                              _:
                                                                  *const fvec_t,
                                                              _: *mut fvec_t)
                                             -> ());
                            current_block = 9437375157805982253;
                        }
                        4 => {
                            (*p).buf = new_fvec(bufsize);
                            (*p).p_object =
                                new_aubio_pitchyinfft(samplerate, bufsize) as
                                    *mut libc::c_void;
                            if (*p).p_object.is_null() {
                                current_block = 9082146709813991249;
                            } else {
                                (*p).detect_cb =
                                    Some(aubio_pitch_do_yinfft as
                                             unsafe extern "C" fn(_:
                                                                      *mut aubio_pitch_t,
                                                                  _:
                                                                      *const fvec_t,
                                                                  _:
                                                                      *mut fvec_t)
                                                 -> ());
                                (*p).conf_cb =
                                    ::std::mem::transmute::<Option<unsafe extern "C" fn(_:
                                                                                            *mut aubio_pitchyinfft_t)
                                                                       ->
                                                                           smpl_t>,
                                                            aubio_pitch_get_conf_t>(Some(aubio_pitchyinfft_get_confidence
                                                                                             as
                                                                                             unsafe extern "C" fn(_:
                                                                                                                      *mut aubio_pitchyinfft_t)
                                                                                                 ->
                                                                                                     smpl_t));
                                aubio_pitchyinfft_set_tolerance((*p).p_object
                                                                    as
                                                                    *mut aubio_pitchyinfft_t,
                                                                0.85f64 as
                                                                    smpl_t);
                                current_block = 9437375157805982253;
                            }
                        }
                        5 => {
                            (*p).buf = new_fvec(bufsize);
                            (*p).p_object =
                                new_aubio_pitchyinfast(bufsize) as
                                    *mut libc::c_void;
                            if (*p).p_object.is_null() {
                                current_block = 9082146709813991249;
                            } else {
                                (*p).detect_cb =
                                    Some(aubio_pitch_do_yinfast as
                                             unsafe extern "C" fn(_:
                                                                      *mut aubio_pitch_t,
                                                                  _:
                                                                      *const fvec_t,
                                                                  _:
                                                                      *mut fvec_t)
                                                 -> ());
                                (*p).conf_cb =
                                    ::std::mem::transmute::<Option<unsafe extern "C" fn(_:
                                                                                            *mut aubio_pitchyinfast_t)
                                                                       ->
                                                                           smpl_t>,
                                                            aubio_pitch_get_conf_t>(Some(aubio_pitchyinfast_get_confidence
                                                                                             as
                                                                                             unsafe extern "C" fn(_:
                                                                                                                      *mut aubio_pitchyinfast_t)
                                                                                                 ->
                                                                                                     smpl_t));
                                aubio_pitchyinfast_set_tolerance((*p).p_object
                                                                     as
                                                                     *mut aubio_pitchyinfast_t,
                                                                 0.15f64 as
                                                                     smpl_t);
                                current_block = 9437375157805982253;
                            }
                        }
                        6 => {
                            (*p).buf = new_fvec(bufsize);
                            (*p).p_object =
                                new_aubio_pitchspecacf(bufsize) as
                                    *mut libc::c_void;
                            if (*p).p_object.is_null() {
                                current_block = 9082146709813991249;
                            } else {
                                (*p).detect_cb =
                                    Some(aubio_pitch_do_specacf as
                                             unsafe extern "C" fn(_:
                                                                      *mut aubio_pitch_t,
                                                                  _:
                                                                      *const fvec_t,
                                                                  _:
                                                                      *mut fvec_t)
                                                 -> ());
                                (*p).conf_cb =
                                    ::std::mem::transmute::<Option<unsafe extern "C" fn(_:
                                                                                            *const aubio_pitchspecacf_t)
                                                                       ->
                                                                           smpl_t>,
                                                            aubio_pitch_get_conf_t>(Some(aubio_pitchspecacf_get_tolerance
                                                                                             as
                                                                                             unsafe extern "C" fn(_:
                                                                                                                      *const aubio_pitchspecacf_t)
                                                                                                 ->
                                                                                                     smpl_t));
                                aubio_pitchspecacf_set_tolerance((*p).p_object
                                                                     as
                                                                     *mut aubio_pitchspecacf_t,
                                                                 0.85f64 as
                                                                     smpl_t);
                                current_block = 9437375157805982253;
                            }
                        }
                        _ => { current_block = 9437375157805982253; }
                    }
                    match current_block {
                        9082146709813991249 => { }
                        _ => { return p }
                    }
                }
            }
        }
    }
    if !(*p).filtered.is_null() { del_fvec((*p).filtered); }
    if !(*p).buf.is_null() { del_fvec((*p).buf); }
    free(p as *mut libc::c_void);
    return 0 as *mut aubio_pitch_t;
}
/* * deletion of the pitch detection object

  \param o pitch detection object as returned by new_aubio_pitch()

*/
#[no_mangle]
pub unsafe extern "C" fn del_aubio_pitch(mut p: *mut aubio_pitch_t) {
    match (*p).type_0 as libc::c_uint {
        0 => {
            del_fvec((*p).buf);
            del_aubio_pitchyin((*p).p_object as *mut aubio_pitchyin_t);
        }
        1 => {
            del_fvec((*p).filtered);
            del_aubio_pvoc((*p).pv);
            del_cvec((*p).fftgrain);
            del_aubio_filter((*p).filter);
            del_aubio_pitchmcomb((*p).p_object as *mut aubio_pitchmcomb_t);
        }
        2 => {
            del_fvec((*p).buf);
            del_aubio_pitchschmitt((*p).p_object as
                                       *mut aubio_pitchschmitt_t);
        }
        3 => {
            del_fvec((*p).buf);
            del_aubio_pitchfcomb((*p).p_object as *mut aubio_pitchfcomb_t);
        }
        4 => {
            del_fvec((*p).buf);
            del_aubio_pitchyinfft((*p).p_object as *mut aubio_pitchyinfft_t);
        }
        5 => {
            del_fvec((*p).buf);
            del_aubio_pitchyinfast((*p).p_object as
                                       *mut aubio_pitchyinfast_t);
        }
        6 => {
            del_fvec((*p).buf);
            del_aubio_pitchspecacf((*p).p_object as
                                       *mut aubio_pitchspecacf_t);
        }
        _ => { }
    }
    free(p as *mut libc::c_void);
}
/* adapter to stack ibuf new samples at the end of buf, and trim `buf` to `bufsize` */
#[no_mangle]
pub unsafe extern "C" fn aubio_pitch_slideblock(mut p: *mut aubio_pitch_t,
                                                mut ibuf: *const fvec_t) {
    let mut overlap_size: uint_t =
        (*(*p).buf).length.wrapping_sub((*ibuf).length);
    // !HAVE_MEMCPY_HACKS
    let mut j: uint_t = 0;
    j = 0 as libc::c_int as uint_t;
    while j < overlap_size {
        *(*(*p).buf).data.offset(j as isize) =
            *(*(*p).buf).data.offset(j.wrapping_add((*ibuf).length) as isize);
        j = j.wrapping_add(1)
    }
    j = 0 as libc::c_int as uint_t;
    while j < (*ibuf).length {
        *(*(*p).buf).data.offset(j.wrapping_add(overlap_size) as isize) =
            *(*ibuf).data.offset(j as isize);
        j = j.wrapping_add(1)
    };
}
/* * set the output unit of the pitch detection object

  \param o pitch detection object as returned by new_aubio_pitch()
  \param mode set pitch units for output

  mode can be one of "Hz", "midi", "cent", or "bin". Defaults to "Hz".

  \return 0 if successfull, non-zero otherwise

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_pitch_set_unit(mut p: *mut aubio_pitch_t,
                                              mut pitch_unit: *const char_t)
 -> uint_t {
    let mut err: uint_t = AUBIO_OK as libc::c_int as uint_t;
    let mut pitch_mode: aubio_pitch_mode = aubio_pitchm_freq;
    if strcmp(pitch_unit, b"freq\x00" as *const u8 as *const libc::c_char) ==
           0 as libc::c_int {
        pitch_mode = aubio_pitchm_freq
    } else if strcmp(pitch_unit,
                     b"hertz\x00" as *const u8 as *const libc::c_char) ==
                  0 as libc::c_int {
        pitch_mode = aubio_pitchm_freq
    } else if strcmp(pitch_unit,
                     b"Hertz\x00" as *const u8 as *const libc::c_char) ==
                  0 as libc::c_int {
        pitch_mode = aubio_pitchm_freq
    } else if strcmp(pitch_unit,
                     b"Hz\x00" as *const u8 as *const libc::c_char) ==
                  0 as libc::c_int {
        pitch_mode = aubio_pitchm_freq
    } else if strcmp(pitch_unit,
                     b"f0\x00" as *const u8 as *const libc::c_char) ==
                  0 as libc::c_int {
        pitch_mode = aubio_pitchm_freq
    } else if strcmp(pitch_unit,
                     b"midi\x00" as *const u8 as *const libc::c_char) ==
                  0 as libc::c_int {
        pitch_mode = aubio_pitchm_midi
    } else if strcmp(pitch_unit,
                     b"cent\x00" as *const u8 as *const libc::c_char) ==
                  0 as libc::c_int {
        pitch_mode = aubio_pitchm_cent
    } else if strcmp(pitch_unit,
                     b"bin\x00" as *const u8 as *const libc::c_char) ==
                  0 as libc::c_int {
        pitch_mode = aubio_pitchm_bin
    } else if strcmp(pitch_unit,
                     b"default\x00" as *const u8 as *const libc::c_char) ==
                  0 as libc::c_int {
        pitch_mode = aubio_pitchm_default
    } else {
        aubio_log(AUBIO_LOG_WRN as libc::c_int,
                  b"AUBIO WARNING: pitch: unknown pitch detection unit \xe2\x80\x98%s\xe2\x80\x99, using default\n\x00"
                      as *const u8 as *const libc::c_char, pitch_unit);
        pitch_mode = aubio_pitchm_default;
        err = AUBIO_FAIL as libc::c_int as uint_t
    }
    (*p).mode = pitch_mode;
    match (*p).mode as libc::c_uint {
        0 => {
            (*p).conv_cb =
                Some(freqconvpass as
                         unsafe extern "C" fn(_: smpl_t, _: uint_t, _: uint_t)
                             -> smpl_t)
        }
        1 => {
            (*p).conv_cb =
                Some(freqconvmidi as
                         unsafe extern "C" fn(_: smpl_t, _: uint_t, _: uint_t)
                             -> smpl_t)
        }
        2 => {
            /* bug: not implemented */
            (*p).conv_cb =
                Some(freqconvmidi as
                         unsafe extern "C" fn(_: smpl_t, _: uint_t, _: uint_t)
                             -> smpl_t)
        }
        3 => {
            (*p).conv_cb =
                Some(freqconvbin as
                         unsafe extern "C" fn(_: smpl_t, _: uint_t, _: uint_t)
                             -> smpl_t)
        }
        _ => { }
    }
    return err;
}
/* * change yin or yinfft tolerance threshold

  \param o pitch detection object as returned by new_aubio_pitch()
  \param tol tolerance default is 0.15 for yin and 0.85 for yinfft

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_pitch_set_tolerance(mut p: *mut aubio_pitch_t,
                                                   mut tol: smpl_t)
 -> uint_t {
    match (*p).type_0 as libc::c_uint {
        0 => {
            aubio_pitchyin_set_tolerance((*p).p_object as
                                             *mut aubio_pitchyin_t, tol);
        }
        4 => {
            aubio_pitchyinfft_set_tolerance((*p).p_object as
                                                *mut aubio_pitchyinfft_t,
                                            tol);
        }
        5 => {
            aubio_pitchyinfast_set_tolerance((*p).p_object as
                                                 *mut aubio_pitchyinfast_t,
                                             tol);
        }
        _ => { }
    }
    return AUBIO_OK as libc::c_int as uint_t;
}
/* * get yin or yinfft tolerance threshold

  \param o pitch detection object as returned by new_aubio_pitch()
  \return tolerance (default is 0.15 for yin and 0.85 for yinfft)

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_pitch_get_tolerance(mut p: *mut aubio_pitch_t)
 -> smpl_t {
    let mut tolerance: smpl_t = 1.0f64 as smpl_t;
    match (*p).type_0 as libc::c_uint {
        0 => {
            tolerance =
                aubio_pitchyin_get_tolerance((*p).p_object as
                                                 *mut aubio_pitchyin_t)
        }
        4 => {
            tolerance =
                aubio_pitchyinfft_get_tolerance((*p).p_object as
                                                    *mut aubio_pitchyinfft_t)
        }
        5 => {
            tolerance =
                aubio_pitchyinfast_get_tolerance((*p).p_object as
                                                     *mut aubio_pitchyinfast_t)
        }
        _ => { }
    }
    return tolerance;
}
/* * set the silence threshold of the pitch detection object

  \param o pitch detection object as returned by new_aubio_pitch()
  \param silence level threshold under which pitch should be ignored, in dB

  \return 0 if successfull, non-zero otherwise

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_pitch_set_silence(mut p: *mut aubio_pitch_t,
                                                 mut silence: smpl_t)
 -> uint_t {
    if silence <= 0 as libc::c_int as libc::c_float &&
           silence >= -(200 as libc::c_int) as libc::c_float {
        (*p).silence = silence;
        return AUBIO_OK as libc::c_int as uint_t
    } else {
        aubio_log(AUBIO_LOG_WRN as libc::c_int,
                  b"AUBIO WARNING: pitch: could not set silence to %.2f\n\x00"
                      as *const u8 as *const libc::c_char,
                  silence as libc::c_double);
        return AUBIO_FAIL as libc::c_int as uint_t
    };
}
/* * set the silence threshold of the pitch detection object

  \param o pitch detection object as returned by ::new_aubio_pitch()

  \return level threshold under which pitch should be ignored, in dB

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_pitch_get_silence(mut p: *mut aubio_pitch_t)
 -> smpl_t {
    return (*p).silence;
}
/* * execute pitch detection on an input signal frame

  \param o pitch detection object as returned by new_aubio_pitch()
  \param in input signal of size [hop_size]
  \param out output pitch candidates of size [1]

*/
/* do method, calling the detection callback, then the conversion callback */
#[no_mangle]
pub unsafe extern "C" fn aubio_pitch_do(mut p: *mut aubio_pitch_t,
                                        mut ibuf: *const fvec_t,
                                        mut obuf: *mut fvec_t) {
    (*p).detect_cb.expect("non-null function pointer")(p, ibuf, obuf);
    if aubio_silence_detection(ibuf, (*p).silence) ==
           1 as libc::c_int as libc::c_uint {
        *(*obuf).data.offset(0 as libc::c_int as isize) = 0.0f64 as smpl_t
    }
    *(*obuf).data.offset(0 as libc::c_int as isize) =
        (*p).conv_cb.expect("non-null function pointer")(*(*obuf).data.offset(0
                                                                                  as
                                                                                  libc::c_int
                                                                                  as
                                                                                  isize),
                                                         (*p).samplerate,
                                                         (*p).bufsize);
}
/* *< silence threshold */
/* callback functions for pitch detection */
/* do method for each algorithm */
unsafe extern "C" fn aubio_pitch_do_mcomb(mut p: *mut aubio_pitch_t,
                                          mut ibuf: *const fvec_t,
                                          mut obuf: *mut fvec_t) {
    aubio_filter_do_outplace((*p).filter, ibuf, (*p).filtered);
    aubio_pvoc_do((*p).pv, ibuf, (*p).fftgrain);
    aubio_pitchmcomb_do((*p).p_object as *mut aubio_pitchmcomb_t,
                        (*p).fftgrain, obuf);
    *(*obuf).data.offset(0 as libc::c_int as isize) =
        aubio_bintofreq(*(*obuf).data.offset(0 as libc::c_int as isize),
                        (*p).samplerate as smpl_t, (*p).bufsize as smpl_t);
}
unsafe extern "C" fn aubio_pitch_do_yin(mut p: *mut aubio_pitch_t,
                                        mut ibuf: *const fvec_t,
                                        mut obuf: *mut fvec_t) {
    let mut pitch: smpl_t = 0.0f64 as smpl_t;
    aubio_pitch_slideblock(p, ibuf);
    aubio_pitchyin_do((*p).p_object as *mut aubio_pitchyin_t, (*p).buf, obuf);
    pitch = *(*obuf).data.offset(0 as libc::c_int as isize);
    if pitch > 0 as libc::c_int as libc::c_float {
        pitch =
            ((*p).samplerate as libc::c_double /
                 (pitch as libc::c_double + 0.0f64)) as smpl_t
    } else { pitch = 0.0f64 as smpl_t }
    *(*obuf).data.offset(0 as libc::c_int as isize) = pitch;
}
unsafe extern "C" fn aubio_pitch_do_yinfft(mut p: *mut aubio_pitch_t,
                                           mut ibuf: *const fvec_t,
                                           mut obuf: *mut fvec_t) {
    let mut pitch: smpl_t = 0.0f64 as smpl_t;
    aubio_pitch_slideblock(p, ibuf);
    aubio_pitchyinfft_do((*p).p_object as *mut aubio_pitchyinfft_t, (*p).buf,
                         obuf);
    pitch = *(*obuf).data.offset(0 as libc::c_int as isize);
    if pitch > 0 as libc::c_int as libc::c_float {
        pitch =
            ((*p).samplerate as libc::c_double /
                 (pitch as libc::c_double + 0.0f64)) as smpl_t
    } else { pitch = 0.0f64 as smpl_t }
    *(*obuf).data.offset(0 as libc::c_int as isize) = pitch;
}
unsafe extern "C" fn aubio_pitch_do_yinfast(mut p: *mut aubio_pitch_t,
                                            mut ibuf: *const fvec_t,
                                            mut obuf: *mut fvec_t) {
    let mut pitch: smpl_t = 0.0f64 as smpl_t;
    aubio_pitch_slideblock(p, ibuf);
    aubio_pitchyinfast_do((*p).p_object as *mut aubio_pitchyinfast_t,
                          (*p).buf, obuf);
    pitch = *(*obuf).data.offset(0 as libc::c_int as isize);
    if pitch > 0 as libc::c_int as libc::c_float {
        pitch =
            ((*p).samplerate as libc::c_double /
                 (pitch as libc::c_double + 0.0f64)) as smpl_t
    } else { pitch = 0.0f64 as smpl_t }
    *(*obuf).data.offset(0 as libc::c_int as isize) = pitch;
}
unsafe extern "C" fn aubio_pitch_do_specacf(mut p: *mut aubio_pitch_t,
                                            mut ibuf: *const fvec_t,
                                            mut out: *mut fvec_t) {
    let mut pitch: smpl_t = 0.0f64 as smpl_t;
    let mut period: smpl_t = 0.;
    aubio_pitch_slideblock(p, ibuf);
    aubio_pitchspecacf_do((*p).p_object as *mut aubio_pitchspecacf_t,
                          (*p).buf, out);
    //out->data[0] = aubio_bintofreq (out->data[0], p->samplerate, p->bufsize);
    period = *(*out).data.offset(0 as libc::c_int as isize);
    if period > 0 as libc::c_int as libc::c_float {
        pitch = (*p).samplerate as libc::c_float / period
    } else { pitch = 0.0f64 as smpl_t }
    *(*out).data.offset(0 as libc::c_int as isize) = pitch;
}
unsafe extern "C" fn aubio_pitch_do_fcomb(mut p: *mut aubio_pitch_t,
                                          mut ibuf: *const fvec_t,
                                          mut out: *mut fvec_t) {
    aubio_pitch_slideblock(p, ibuf);
    aubio_pitchfcomb_do((*p).p_object as *mut aubio_pitchfcomb_t, (*p).buf,
                        out);
    *(*out).data.offset(0 as libc::c_int as isize) =
        aubio_bintofreq(*(*out).data.offset(0 as libc::c_int as isize),
                        (*p).samplerate as smpl_t, (*p).bufsize as smpl_t);
}
unsafe extern "C" fn aubio_pitch_do_schmitt(mut p: *mut aubio_pitch_t,
                                            mut ibuf: *const fvec_t,
                                            mut out: *mut fvec_t) {
    let mut period: smpl_t = 0.;
    let mut pitch: smpl_t = 0.0f64 as smpl_t;
    aubio_pitch_slideblock(p, ibuf);
    aubio_pitchschmitt_do((*p).p_object as *mut aubio_pitchschmitt_t,
                          (*p).buf, out);
    period = *(*out).data.offset(0 as libc::c_int as isize);
    if period > 0 as libc::c_int as libc::c_float {
        pitch = (*p).samplerate as libc::c_float / period
    } else { pitch = 0.0f64 as smpl_t }
    *(*out).data.offset(0 as libc::c_int as isize) = pitch;
}
/* internal functions for frequency conversion */
/* conversion callbacks */
unsafe extern "C" fn freqconvbin(mut f: smpl_t, mut samplerate: uint_t,
                                 mut bufsize: uint_t) -> smpl_t {
    return aubio_freqtobin(f, samplerate as smpl_t, bufsize as smpl_t);
}
unsafe extern "C" fn freqconvmidi(mut f: smpl_t, mut samplerate: uint_t,
                                  mut bufsize: uint_t) -> smpl_t {
    return aubio_freqtomidi(f);
}
unsafe extern "C" fn freqconvpass(mut f: smpl_t, mut samplerate: uint_t,
                                  mut bufsize: uint_t) -> smpl_t {
    return f;
}
/* * get the current confidence

  \param o pitch detection object as returned by new_aubio_pitch()

  \return the current confidence of the pitch algorithm

*/
/* confidence callbacks */
#[no_mangle]
pub unsafe extern "C" fn aubio_pitch_get_confidence(mut p: *mut aubio_pitch_t)
 -> smpl_t {
    if (*p).conf_cb.is_some() {
        return (*p).conf_cb.expect("non-null function pointer")((*p).p_object)
    }
    return 0.0f64 as smpl_t;
}
