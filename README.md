# GNSS-SDRLIB-PVT
Modified version of GNSS-SDRLIB

(Note: This readme document is in work).

Purpose. This is a significantly-modifed version of GNSS-SDRLIB. Instead of requiring a TCP link to RTKLIB for positioning, it has its own least squares solver built-in. 

Limitations. This is a highly stripped-down version of GNSS-SDRLIB, in an attempt to create a more truly real-time and robust GPS SDR. This version is GPS L1 only, Linux only (not Windows), and only currently works with RTLSDR or bladeRF front ends.

Videos. The following three videos are available on Youtube:

Part 1. Overview of Features and Modifications. Link: https://youtu.be/VLlu77O3vec

Part 2. Hardware Used, and Software Builds. Link: https://youtu.be/R3IFhFoblQU 

Part 3. Demo of Playback and Real-Time Modes. Link: https://youtu.be/8pV1hX6IkFk 

Acknowledgements. This software is a modified version of the original 2014 GNSS-SDRLIB code (https://github.com/taroz/GNSS-SDRLIB). Please refer to the videos listed above for a description of features and the modifications made. The author is very grateful for Taro Suzuki's original development of GNSS-SDRLIB, as it is a terrific learning tool, is thoughtfully written, with clear variables and code comments. 

License. This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Keywords. RTLSDR, bladeRF, SDR, software defined radio, GPS, GNSS, L1, GNSS-SDRLIB, GNSS-SDRLIB-PVT, PVT, C, SDL2, NML.

