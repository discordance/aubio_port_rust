use crate::transpiled::awhitening::aubio_spectral_whitening_t;
use crate::transpiled::peakpicker::aubio_peakpicker_t;
use crate::transpiled::phasevoc::aubio_pvoc_t;
use crate::transpiled::specdesc::aubio_specdesc_t;

extern "C" {
    // pub type _aubio_specdesc_t;
    // pub type _aubio_spectral_whitening_t;
    // pub type _aubio_peakpicker_t;
    #[no_mangle]
    fn floorf(_: f32) -> f32;
    #[no_mangle]
    fn free(_: *mut core::ffi::c_void);
    #[no_mangle]
    fn calloc(_: u64, _: u64) -> *mut core::ffi::c_void;
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
    /* * check if buffer level in dB SPL is under a given threshold

      \param v vector to get level from
      \param threshold threshold in dB SPL

      \return 0 if level is under the given threshold, 1 otherwise

    */
    #[no_mangle]
    fn aubio_silence_detection(v: *const fvec_t, threshold: smpl_t) -> uint_t;
    #[no_mangle]
    fn aubio_log(level: sint_t, fmt: *const char_t, _: ...) -> uint_t;
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
    /* * take logarithmic magnitude

      \param s input cvec to compress
      \param lambda value to use for normalisation

      \f$ S_k = log( \lambda * S_k + 1 ) \f$

    */
    #[no_mangle]
    fn cvec_logmag(s: *mut cvec_t, lambda: smpl_t);
    /* * execute spectral description function on a spectral frame

      Generic function to compute spectral description.

      \param o spectral description object as returned by new_aubio_specdesc()
      \param fftgrain input signal spectrum as computed by aubio_pvoc_do
      \param desc output vector (one sample long, to send to the peak picking)

    */
    #[no_mangle]
    fn aubio_specdesc_do(o: *mut aubio_specdesc_t, fftgrain: *const cvec_t, desc: *mut fvec_t);
    /* * creation of a spectral description object

      \param method spectral description method
      \param buf_size length of the input spectrum frame

      The parameter \p method is a string that can be any of:

        - onset novelty functions: `complex`, `energy`, `hfc`, `kl`, `mkl`,
        `phase`, `specdiff`, `specflux`, `wphase`,

        - spectral descriptors: `centroid`, `decrease`, `kurtosis`, `rolloff`,
        `skewness`, `slope`, `spread`.

    */
    #[no_mangle]
    fn new_aubio_specdesc(method: *const char_t, buf_size: uint_t) -> *mut aubio_specdesc_t;
    /* * deletion of a spectral descriptor

      \param o spectral descriptor object as returned by new_aubio_specdesc()

    */
    #[no_mangle]
    fn del_aubio_specdesc(o: *mut aubio_specdesc_t);
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
    fn aubio_pvoc_do(pv: *mut aubio_pvoc_t, in_0: *const fvec_t, fftgrain: *mut cvec_t);
    /* * execute spectral adaptive whitening, in-place

      \param o spectral whitening object as returned by new_aubio_spectral_whitening()
      \param fftgrain input signal spectrum as computed by aubio_pvoc_do() or aubio_fft_do()

    */
    #[no_mangle]
    fn aubio_spectral_whitening_do(o: *mut aubio_spectral_whitening_t, fftgrain: *mut cvec_t);
    /* * creation of a spectral whitening object

      \param buf_size window size of input grains
      \param hop_size number of samples between two consecutive input grains
      \param samplerate sampling rate of the input signal

    */
    #[no_mangle]
    fn new_aubio_spectral_whitening(
        buf_size: uint_t,
        hop_size: uint_t,
        samplerate: uint_t,
    ) -> *mut aubio_spectral_whitening_t;
    /* * set relaxation time for spectral whitening

    \param o spectral whitening object as returned by new_aubio_spectral_whitening()
    \param relax_time relaxation time in seconds between 20 and 500, defaults 250

    */
    #[no_mangle]
    fn aubio_spectral_whitening_set_relax_time(
        o: *mut aubio_spectral_whitening_t,
        relax_time: smpl_t,
    ) -> uint_t;
    /* * set floor for spectral whitening

    \param o spectral whitening object as returned by new_aubio_spectral_whitening()
    \param floor value (typically between 1.e-6 and .2, defaults to 1.e-4)

    */
    #[no_mangle]
    fn aubio_spectral_whitening_set_floor(
        o: *mut aubio_spectral_whitening_t,
        floor: smpl_t,
    ) -> uint_t;
    /* * deletion of a spectral whitening

      \param o spectral whitening object as returned by new_aubio_spectral_whitening()

    */
    #[no_mangle]
    fn del_aubio_spectral_whitening(o: *mut aubio_spectral_whitening_t);
    /* * peak-picker creation function */
    #[no_mangle]
    fn new_aubio_peakpicker() -> *mut aubio_peakpicker_t;
    /* * real time peak picking function */
    #[no_mangle]
    fn aubio_peakpicker_do(p: *mut aubio_peakpicker_t, in_0: *mut fvec_t, out: *mut fvec_t);
    /* * destroy peak picker structure */
    #[no_mangle]
    fn del_aubio_peakpicker(p: *mut aubio_peakpicker_t);
    /* * get current peak value */
    #[no_mangle]
    fn aubio_peakpicker_get_thresholded_input(p: *mut aubio_peakpicker_t) -> *mut fvec_t;
    /* * set peak picking threshold */
    #[no_mangle]
    fn aubio_peakpicker_set_threshold(p: *mut aubio_peakpicker_t, threshold: smpl_t) -> uint_t;
    /* * get peak picking threshold */
    #[no_mangle]
    fn aubio_peakpicker_get_threshold(p: *mut aubio_peakpicker_t) -> smpl_t;
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
pub type C2RustUnnamed = u32;
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
pub type aubio_log_level = u32;
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

  Spectral adaptive whitening

  References:

  D. Stowell and M. D. Plumbley. Adaptive whitening for improved real-time
  audio onset detection. In Proceedings of the International Computer Music
  Conference (ICMC), 2007, Copenhagen, Denmark.

  http://www.eecs.qmul.ac.uk/~markp/2007/StowellPlumbley07-icmc.pdf

  S. Böck,, F. Krebs, and M. Schedl. Evaluating the Online Capabilities of
  Onset Detection Methods. In Proceedings of the 13th International Society for
  Music Information Retrieval Conference (ISMIR), 2012, Porto, Portugal.

  http://ismir2012.ismir.net/event/papers/049_ISMIR_2012.pdf
  http://www.cp.jku.at/research/papers/Boeck_etal_ISMIR_2012.pdf

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

  Peak picking utilities function

  \example onset/test-peakpicker.c

*/
/* * structure to store object state */
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _aubio_onset_t {
    pub pv: *mut aubio_pvoc_t,
    pub od: *mut aubio_specdesc_t,
    pub pp: *mut aubio_peakpicker_t,
    pub fftgrain: *mut cvec_t,
    pub desc: *mut fvec_t,
    pub silence: smpl_t,
    pub minioi: uint_t,
    pub delay: uint_t,
    pub samplerate: uint_t,
    pub hop_size: uint_t,
    pub total_frames: uint_t,
    pub last_onset: uint_t,
    pub apply_compression: uint_t,
    pub lambda_compression: smpl_t,
    pub apply_awhitening: uint_t,
    pub spectral_whitening: *mut aubio_spectral_whitening_t,
}
/*
  Copyright (C) 2006-2013 Paul Brossier <piem@aubio.org>

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

  Onset detection object

  The following routines compute the onset detection function and detect peaks
  in these functions. When onsets are found above a given silence threshold,
  and after a minimum inter-onset interval, the output vector returned by
  aubio_onset_do() is filled with `1`. Otherwise, the output vector remains
  `0`.

  The peak-picking threshold, the silence threshold, and the minimum
  inter-onset interval can be adjusted during the execution of the
  aubio_onset_do routine using the corresponding functions.

  \example onset/test-onset.c
  \example examples/aubioonset.c
  \example examples/aubionotes.c

*/
/* * onset detection object */
pub type aubio_onset_t = _aubio_onset_t;
/* * execute onset detection

  \param o onset detection object as returned by new_aubio_onset()
  \param input new audio vector of length hop_size
  \param onset output vector of length 1, containing 0 if no onset was found,
  and a value equal or greater than 1 otherwise

  When no onset was detected, the first element of the output vector `onset`
  is set to 0.

  When an onset is found, the first element of the output vector `onset` is set
  to `offset = 1 + a` where `a` is a number in the range`[0, 1]`.

  The final onset detection time, in samples, can be obtained with
  aubio_onset_get_last(). It can also be derived from `offset` as
  follows:

  \code
    t = total_frames + offset * hop_size - delay
  \endcode

  where `total_frames` is the total number of frames processed so far, and
  `delay` is the current delay of the onset object, as returned by
  aubio_onset_get_delay().

*/
/* execute onset detection function on iput buffer */
#[no_mangle]
pub unsafe extern "C" fn aubio_onset_do(
    mut o: *mut aubio_onset_t,
    mut input: *const fvec_t,
    mut onset: *mut fvec_t,
) {
    let mut isonset: smpl_t = 0 as i32 as smpl_t;
    aubio_pvoc_do((*o).pv, input, (*o).fftgrain);
    /*
    if (apply_filtering) {
    }
    */
    if (*o).apply_awhitening != 0 {
        aubio_spectral_whitening_do((*o).spectral_whitening, (*o).fftgrain);
    }
    if (*o).apply_compression != 0 {
        cvec_logmag((*o).fftgrain, (*o).lambda_compression);
    }
    aubio_specdesc_do((*o).od, (*o).fftgrain, (*o).desc);
    aubio_peakpicker_do((*o).pp, (*o).desc, onset);
    isonset = *(*onset).data.offset(0 as i32 as isize);
    if isonset as f64 > 0.0f64 {
        if aubio_silence_detection(input, (*o).silence) == 1 as i32 as u32 {
            //AUBIO_DBG ("silent onset, not marking as onset\n");
            isonset = 0 as i32 as smpl_t
        } else {
            // we have an onset
            let mut new_onset: uint_t = (*o).total_frames.wrapping_add(floorf(
                ((isonset * (*o).hop_size as f32) as f64 + 0.5f64) as f32,
            ) as uint_t);
            // check if last onset time was more than minioi ago
            if (*o).last_onset.wrapping_add((*o).minioi) < new_onset {
                // start of file: make sure (new_onset - delay) >= 0
                if (*o).last_onset > 0 as i32 as u32 && (*o).delay > new_onset {
                    isonset = 0 as i32 as smpl_t
                } else {
                    //AUBIO_DBG ("accepted detection, marking as onset\n");
                    (*o).last_onset = if (*o).delay > new_onset {
                        (*o).delay
                    } else {
                        new_onset
                    }
                }
            } else {
                //AUBIO_DBG ("doubled onset, not marking as onset\n");
                isonset = 0 as i32 as smpl_t
            }
        }
    } else if (*o).total_frames <= (*o).delay {
        // we are at the beginning of the file
        // and we don't find silence
        if aubio_silence_detection(input, (*o).silence) == 0 as i32 as u32 {
            let mut new_onset_0: uint_t = (*o).total_frames;
            if (*o).total_frames == 0 as i32 as u32
                || (*o).last_onset.wrapping_add((*o).minioi) < new_onset_0
            {
                isonset = (*o).delay.wrapping_div((*o).hop_size) as smpl_t;
                (*o).last_onset = (*o).total_frames.wrapping_add((*o).delay)
            }
        }
    }
    *(*onset).data.offset(0 as i32 as isize) = isonset;
    (*o).total_frames = ((*o).total_frames as u32).wrapping_add((*o).hop_size) as uint_t as uint_t;
}
/* * get the time of the latest onset detected, in samples

  \param o onset detection object as returned by new_aubio_onset()

  \return onset detection timestamps (in samples)

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_onset_get_last(mut o: *const aubio_onset_t) -> uint_t {
    return (*o).last_onset.wrapping_sub((*o).delay);
}
/* * get the time of the latest onset detected, in seconds

  \param o onset detection object as returned by new_aubio_onset()

  \return onset detection timestamps (in seconds)

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_onset_get_last_s(mut o: *const aubio_onset_t) -> smpl_t {
    return aubio_onset_get_last(o) as f32 / (*o).samplerate as smpl_t;
}
/* * get the time of the latest onset detected, in milliseconds

  \param o onset detection object as returned by new_aubio_onset()

  \return onset detection timestamps (in milliseconds)

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_onset_get_last_ms(mut o: *const aubio_onset_t) -> smpl_t {
    return (aubio_onset_get_last_s(o) as f64 * 1000.0f64) as smpl_t;
}
/* * set onset detection adaptive whitening

  \param o onset detection object as returned by new_aubio_onset()
  \param enable 1 to enable, 0 to disable

  \return 0 if successful, 1 otherwise

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_onset_set_awhitening(
    mut o: *mut aubio_onset_t,
    mut enable: uint_t,
) -> uint_t {
    (*o).apply_awhitening = if enable == 1 as i32 as u32 {
        1 as i32
    } else {
        0 as i32
    } as uint_t;
    return AUBIO_OK as i32 as uint_t;
}
/* * get onset detection adaptive whitening

  \param o onset detection object as returned by new_aubio_onset()

  \return 1 if enabled, 0 otherwise

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_onset_get_awhitening(mut o: *mut aubio_onset_t) -> smpl_t {
    return (*o).apply_awhitening as smpl_t;
}
/* * set or disable log compression

 \param o onset detection object as returned by new_aubio_onset()
 \param lambda logarithmic compression factor, 0 to disable

 \return 0 if successful, 1 otherwise

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_onset_set_compression(
    mut o: *mut aubio_onset_t,
    mut lambda: smpl_t,
) -> uint_t {
    if (lambda as f64) < 0.0f64 {
        return AUBIO_FAIL as i32 as uint_t;
    }
    (*o).lambda_compression = lambda;
    (*o).apply_compression = if (*o).lambda_compression as f64 > 0.0f64 {
        1 as i32
    } else {
        0 as i32
    } as uint_t;
    return AUBIO_OK as i32 as uint_t;
}
/* * get onset detection log compression

 \param o onset detection object as returned by new_aubio_onset()

 \returns 0 if disabled, compression factor otherwise

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_onset_get_compression(mut o: *mut aubio_onset_t) -> smpl_t {
    return if (*o).apply_compression != 0 {
        (*o).lambda_compression
    } else {
        0 as i32 as f32
    };
}
/* * set onset detection silence threshold

  \param o onset detection object as returned by new_aubio_onset()
  \param silence new silence detection threshold

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_onset_set_silence(
    mut o: *mut aubio_onset_t,
    mut silence: smpl_t,
) -> uint_t {
    (*o).silence = silence;
    return AUBIO_OK as i32 as uint_t;
}
/* * get onset detection silence threshold

  \param o onset detection object as returned by new_aubio_onset()

  \return current silence threshold

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_onset_get_silence(mut o: *const aubio_onset_t) -> smpl_t {
    return (*o).silence;
}
/* * set onset detection peak picking threshold

  \param o onset detection object as returned by new_aubio_onset()
  \param threshold new peak-picking threshold

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_onset_set_threshold(
    mut o: *mut aubio_onset_t,
    mut threshold: smpl_t,
) -> uint_t {
    aubio_peakpicker_set_threshold((*o).pp, threshold);
    return AUBIO_OK as i32 as uint_t;
}
/* * get onset peak picking threshold

  \param o onset detection object as returned by new_aubio_onset()
  \return current onset detection threshold

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_onset_get_threshold(mut o: *const aubio_onset_t) -> smpl_t {
    return aubio_peakpicker_get_threshold((*o).pp);
}
/* * set minimum inter onset interval in samples

  \param o onset detection object as returned by new_aubio_onset()
  \param minioi minimum interval between two consecutive onsets (in
  samples)

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_onset_set_minioi(
    mut o: *mut aubio_onset_t,
    mut minioi: uint_t,
) -> uint_t {
    (*o).minioi = minioi;
    return AUBIO_OK as i32 as uint_t;
}
/* * get minimum inter onset interval in samples

  \param o onset detection object as returned by new_aubio_onset()
  \return minimum interval between two consecutive onsets (in
  samples)

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_onset_get_minioi(mut o: *const aubio_onset_t) -> uint_t {
    return (*o).minioi;
}
/* * set minimum inter onset interval in seconds

  \param o onset detection object as returned by new_aubio_onset()
  \param minioi minimum interval between two consecutive onsets (in
  seconds)

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_onset_set_minioi_s(
    mut o: *mut aubio_onset_t,
    mut minioi: smpl_t,
) -> uint_t {
    return aubio_onset_set_minioi(
        o,
        floorf(((minioi * (*o).samplerate as f32) as f64 + 0.5f64) as f32) as uint_t,
    );
}
/* * get minimum inter onset interval in seconds

  \param o onset detection object as returned by new_aubio_onset()
  \return minimum interval between two consecutive onsets (in
  seconds)

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_onset_get_minioi_s(mut o: *const aubio_onset_t) -> smpl_t {
    return aubio_onset_get_minioi(o) as f32 / (*o).samplerate as smpl_t;
}
/* * set minimum inter onset interval in milliseconds

  \param o onset detection object as returned by new_aubio_onset()
  \param minioi minimum interval between two consecutive onsets (in
  milliseconds)

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_onset_set_minioi_ms(
    mut o: *mut aubio_onset_t,
    mut minioi: smpl_t,
) -> uint_t {
    return aubio_onset_set_minioi_s(o, (minioi as f64 / 1000.0f64) as smpl_t);
}
/* * get minimum inter onset interval in milliseconds

  \param o onset detection object as returned by new_aubio_onset()
  \return minimum interval between two consecutive onsets (in
  milliseconds)

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_onset_get_minioi_ms(mut o: *const aubio_onset_t) -> smpl_t {
    return (aubio_onset_get_minioi_s(o) as f64 * 1000.0f64) as smpl_t;
}
/* * set delay in samples

  \param o onset detection object as returned by new_aubio_onset()
  \param delay constant system delay to take back from detection time
  (in samples)

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_onset_set_delay(
    mut o: *mut aubio_onset_t,
    mut delay: uint_t,
) -> uint_t {
    (*o).delay = delay;
    return AUBIO_OK as i32 as uint_t;
}
/* * get delay in samples

  \param o onset detection object as returned by new_aubio_onset()
  \return constant system delay to take back from detection time
  (in samples)

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_onset_get_delay(mut o: *const aubio_onset_t) -> uint_t {
    return (*o).delay;
}
/* * set delay in seconds

  \param o onset detection object as returned by new_aubio_onset()
  \param delay constant system delay to take back from detection time
  (in seconds)

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_onset_set_delay_s(
    mut o: *mut aubio_onset_t,
    mut delay: smpl_t,
) -> uint_t {
    return aubio_onset_set_delay(o, (delay * (*o).samplerate as f32) as uint_t);
}
/* * get delay in seconds

  \param o onset detection object as returned by new_aubio_onset()
  \return constant system delay to take back from detection time
  (in seconds)

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_onset_get_delay_s(mut o: *const aubio_onset_t) -> smpl_t {
    return aubio_onset_get_delay(o) as f32 / (*o).samplerate as smpl_t;
}
/* * set delay in milliseconds

  \param o onset detection object as returned by new_aubio_onset()
  \param delay constant system delay to take back from detection time
  (in milliseconds)

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_onset_set_delay_ms(
    mut o: *mut aubio_onset_t,
    mut delay: smpl_t,
) -> uint_t {
    return aubio_onset_set_delay_s(o, (delay as f64 / 1000.0f64) as smpl_t);
}
/* * get delay in milliseconds

  \param o onset detection object as returned by new_aubio_onset()
  \return constant system delay to take back from detection time
  (in milliseconds)

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_onset_get_delay_ms(mut o: *const aubio_onset_t) -> smpl_t {
    return (aubio_onset_get_delay_s(o) as f64 * 1000.0f64) as smpl_t;
}
/* * get onset detection function

  \param o onset detection object as returned by new_aubio_onset()
  \return the current value of the descriptor

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_onset_get_descriptor(mut o: *const aubio_onset_t) -> smpl_t {
    return *(*(*o).desc).data.offset(0 as i32 as isize);
}
/* * get thresholded onset detection function

  \param o onset detection object as returned by new_aubio_onset()
  \return the value of the thresholded descriptor

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_onset_get_thresholded_descriptor(
    mut o: *const aubio_onset_t,
) -> smpl_t {
    let mut thresholded: *mut fvec_t = aubio_peakpicker_get_thresholded_input((*o).pp);
    return *(*thresholded).data.offset(0 as i32 as isize);
}
/* * create onset detection object

  \param method onset detection type as specified in specdesc.h
  \param buf_size buffer size for phase vocoder
  \param hop_size hop size for phase vocoder
  \param samplerate sampling rate of the input signal

  \return newly created ::aubio_onset_t

*/
/* Allocate memory for an onset detection */
#[no_mangle]
pub unsafe extern "C" fn new_aubio_onset(
    mut onset_mode: *const char_t,
    mut buf_size: uint_t,
    mut hop_size: uint_t,
    mut samplerate: uint_t,
) -> *mut aubio_onset_t {
    let mut o: *mut aubio_onset_t = calloc(
        ::std::mem::size_of::<aubio_onset_t>() as u64,
        1 as i32 as u64,
    ) as *mut aubio_onset_t;
    /* check parameters are valid */
    if (hop_size as sint_t) < 1 as i32 {
        aubio_log(
            AUBIO_LOG_ERR as i32,
            b"AUBIO ERROR: onset: got hop_size %d, but can not be < 1\n\x00" as *const u8
                as *const i8,
            hop_size,
        );
    } else if (buf_size as sint_t) < 2 as i32 {
        aubio_log(
            AUBIO_LOG_ERR as i32,
            b"AUBIO ERROR: onset: got buffer_size %d, but can not be < 2\n\x00" as *const u8
                as *const i8,
            buf_size,
        );
    } else if buf_size < hop_size {
        aubio_log(
            AUBIO_LOG_ERR as i32,
            b"AUBIO ERROR: onset: hop size (%d) is larger than win size (%d)\n\x00" as *const u8
                as *const i8,
            hop_size,
            buf_size,
        );
    } else if (samplerate as sint_t) < 1 as i32 {
        aubio_log(
            AUBIO_LOG_ERR as i32,
            b"AUBIO ERROR: onset: samplerate (%d) can not be < 1\n\x00" as *const u8 as *const i8,
            samplerate,
        );
    } else {
        /* store creation parameters */
        (*o).samplerate = samplerate;
        (*o).hop_size = hop_size;
        /* allocate memory */
        (*o).pv = new_aubio_pvoc(buf_size, (*o).hop_size);
        (*o).pp = new_aubio_peakpicker();
        (*o).od = new_aubio_specdesc(onset_mode, buf_size);
        (*o).fftgrain = new_cvec(buf_size);
        (*o).desc = new_fvec(1 as i32 as uint_t);
        (*o).spectral_whitening = new_aubio_spectral_whitening(buf_size, hop_size, samplerate);
        if !((*o).pv.is_null()
            || (*o).pp.is_null()
            || (*o).od.is_null()
            || (*o).fftgrain.is_null()
            || (*o).desc.is_null()
            || (*o).spectral_whitening.is_null())
        {
            /* initialize internal variables */
            aubio_onset_set_default_parameters(o, onset_mode);
            aubio_onset_reset(o);
            return o;
        }
    }
    del_aubio_onset(o);
    return 0 as *mut aubio_onset_t;
}
/* * reset onset detection

 \param o onset detection object as returned by new_aubio_onset()

 Reset current time and last onset to 0.

 This function is called at the end of new_aubio_onset().

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_onset_reset(mut o: *mut aubio_onset_t) {
    (*o).last_onset = 0 as i32 as uint_t;
    (*o).total_frames = 0 as i32 as uint_t;
}
/* * set default parameters

 \param o onset detection object as returned by new_aubio_onset()
 \param onset_mode detection mode to adjust

 This function is called at the end of new_aubio_onset().

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_onset_set_default_parameters(
    mut o: *mut aubio_onset_t,
    mut onset_mode: *const char_t,
) -> uint_t {
    let mut ret: uint_t = AUBIO_OK as i32 as uint_t;
    /* set some default parameter */
    aubio_onset_set_threshold(o, 0.3f64 as smpl_t);
    aubio_onset_set_delay(o, (4.3f64 * (*o).hop_size as f64) as uint_t);
    aubio_onset_set_minioi_ms(o, 50.0f64 as smpl_t);
    aubio_onset_set_silence(o, -70.0f64 as smpl_t);
    // disable spectral whitening
    aubio_onset_set_awhitening(o, 0 as i32 as uint_t);
    // disable logarithmic magnitude
    aubio_onset_set_compression(o, 0.0f64 as smpl_t);
    /* method specific optimisations */
    if !(strcmp(onset_mode, b"energy\x00" as *const u8 as *const i8) == 0 as i32) {
        if strcmp(onset_mode, b"hfc\x00" as *const u8 as *const i8) == 0 as i32
            || strcmp(onset_mode, b"default\x00" as *const u8 as *const i8) == 0 as i32
        {
            aubio_onset_set_threshold(o, 0.058f64 as smpl_t);
            aubio_onset_set_compression(o, 1.0f64 as smpl_t);
        } else if strcmp(onset_mode, b"complexdomain\x00" as *const u8 as *const i8) == 0 as i32
            || strcmp(onset_mode, b"complex\x00" as *const u8 as *const i8) == 0 as i32
        {
            aubio_onset_set_delay(o, (4.6f64 * (*o).hop_size as f64) as uint_t);
            aubio_onset_set_threshold(o, 0.15f64 as smpl_t);
            aubio_onset_set_awhitening(o, 1 as i32 as uint_t);
            aubio_onset_set_compression(o, 1.0f64 as smpl_t);
        } else if strcmp(onset_mode, b"phase\x00" as *const u8 as *const i8) == 0 as i32 {
            (*o).apply_compression = 0 as i32 as uint_t;
            aubio_onset_set_awhitening(o, 0 as i32 as uint_t);
        } else if !(strcmp(onset_mode, b"wphase\x00" as *const u8 as *const i8) == 0 as i32) {
            if strcmp(onset_mode, b"mkl\x00" as *const u8 as *const i8) == 0 as i32 {
                aubio_onset_set_threshold(o, 0.05f64 as smpl_t);
                aubio_onset_set_awhitening(o, 1 as i32 as uint_t);
                aubio_onset_set_compression(o, 0.02f64 as smpl_t);
            } else if strcmp(onset_mode, b"kl\x00" as *const u8 as *const i8) == 0 as i32 {
                aubio_onset_set_threshold(o, 0.35f64 as smpl_t);
                aubio_onset_set_awhitening(o, 1 as i32 as uint_t);
                aubio_onset_set_compression(o, 0.02f64 as smpl_t);
            } else if strcmp(onset_mode, b"specflux\x00" as *const u8 as *const i8) == 0 as i32 {
                aubio_onset_set_threshold(o, 0.18f64 as smpl_t);
                aubio_onset_set_awhitening(o, 1 as i32 as uint_t);
                aubio_spectral_whitening_set_relax_time(
                    (*o).spectral_whitening,
                    100 as i32 as smpl_t,
                );
                aubio_spectral_whitening_set_floor((*o).spectral_whitening, 1.0f64 as smpl_t);
                aubio_onset_set_compression(o, 10.0f64 as smpl_t);
            } else if !(strcmp(onset_mode, b"specdiff\x00" as *const u8 as *const i8) == 0 as i32) {
                if strcmp(onset_mode, b"old_default\x00" as *const u8 as *const i8) == 0 as i32 {
                    // used to reproduce results obtained with the previous version
                    aubio_onset_set_threshold(o, 0.3f64 as smpl_t);
                    aubio_onset_set_minioi_ms(o, 20.0f64 as smpl_t);
                    aubio_onset_set_compression(o, 0.0f64 as smpl_t);
                } else {
                    aubio_log(AUBIO_LOG_WRN as i32,
                              b"AUBIO WARNING: onset: unknown spectral descriptor type %s, using default parameters.\n\x00"
                                  as *const u8 as *const i8,
                              onset_mode);
                    ret = AUBIO_FAIL as i32 as uint_t
                }
            }
        }
    }
    return ret;
}
/* * delete onset detection object

  \param o onset detection object to delete

*/
#[no_mangle]
pub unsafe extern "C" fn del_aubio_onset(mut o: *mut aubio_onset_t) {
    if !(*o).spectral_whitening.is_null() {
        del_aubio_spectral_whitening((*o).spectral_whitening);
    }
    if !(*o).od.is_null() {
        del_aubio_specdesc((*o).od);
    }
    if !(*o).pp.is_null() {
        del_aubio_peakpicker((*o).pp);
    }
    if !(*o).pv.is_null() {
        del_aubio_pvoc((*o).pv);
    }
    if !(*o).desc.is_null() {
        del_fvec((*o).desc);
    }
    if !(*o).fftgrain.is_null() {
        del_cvec((*o).fftgrain);
    }
    free(o as *mut core::ffi::c_void);
}
