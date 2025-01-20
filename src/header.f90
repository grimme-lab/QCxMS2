! This file is part of QCxMS2.
!
! SPDX-Identifier: LGPL-3.0-or-later
!
! QCxMS2 is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! QCxMS2 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with QCxMS2.  If not, see <https://www.gnu.org/licenses/>.

subroutine qcxms2_header

   write (*, '(a)') &
      !< < < < < < < < < < < < < < < < < < > > > > > > > > > > > > > > > > > >!
      "      -----------------------------------------------------------      ", &
      "     |                   =====================                   |     ", &
      "     |                   =       QCxMS2      =                   |     ", &
      "     |                   =====================                   |     ", &
      "     |                      J. Gorges                            |     ", &
      "     |                      S. Grimme                            |     ", &
      "     |          Mulliken Center for Theoretical Chemistry        |     ", &
      "     |                    University of Bonn                     |     ", &
      !   "     |  Version number by <",version," (<",author,"><",date,">)  |     ", &
      "      -----------------------------------------------------------      ", ""
   !< < < < < < < < < < < < < < < < < < > > > > > > > > > > > > > > > > > >!
   call qcxms2_version

end subroutine qcxms2_header

subroutine qcxms2_version

   include 'qcxms2_version.fh'
   write (*, '(3x,"*",*(1x,a))') &
      & "QCxMS2 version", version, "compiled by", author, "on", date
end subroutine qcxms2_version

subroutine disclamer()

   write (*, '(3x,a)') &
      !< < < < < < < < < < < < < < < > > > > > > > > > > > > > > > >!
      "QCxMS2 is free software: you can redistribute it and/or modify it under", &
      "the terms of the GNU Lesser General Public License as published by", &
      "the Free Software Foundation, either version 3 of the License, or", &
      "(at your option) any later version.", &
      "", &
      "QCxMS2 is distributed in the hope that it will be useful,", &
      "but WITHOUT ANY WARRANTY; without even the implied warranty of", &
      "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the", &
      "GNU Lesser General Public License for more details.", &
      ""
end subroutine disclamer

subroutine qcxms2_logo
   write (*, '(a)') &
      "-----------------------------------------------------------------------------------------------------", &
      "|d000Okoc:;;:lxO0O; 'k000Oxoc:;;:ldO0000000000000000Oxllldk00000000kollok00000koc;;;:lxO00koc;;:cokOd|", &
      "|OMWKo'.      .:kNl ;XMW0c'.      .;0MMMMMMMMMMMMMMMNc   .dWMMMMMMXc.  .xMMMKo'.      .d0o,.     .;kk|", &
      "|OWx.    .'.    .c; ;XNd.    .'..  ,0MMMMMMMMMMMMMMMX;    .kWMMMMWo     lWM0'    .''. .k0;  .''.   .;|", &
      "|Ox.   .dKNXk,      ;Xd.   'xKNNKkokOxxxkXMMMKxxxxKM0'     ,KMMMMk.     :NNc   .oXNNKkOWMXdo0NNx.   .|", &
      "|d;   .dWMMMMK,     :O,   .kMMMMMMMWx.   :KWO'   ;KMk.      cNMM0,      ;XX:   .xWMMMMMMMMMMMMMk.   .|", &
      "|c.   ,KMMMMMWl     :d.   ;XMMMMMMMMWx.   :d'   ;KMMx.      .dWNc       ,KWx.   .lONMMMMMMMMMMX:   .;|", &
      "|,    :NMMMMMMx.    ;l.   :NMMMMMMMMMWk.       :KMMMd   .'.  .Od.  ..   '0MNx'    .,oKWMMMMMW0;   .ok|", &
      "|,    :NMMMMMMx.    :l.   :NMMMMMMMMMMWk.     :KMMMWo   .d:   ..   ol   .OMMMXd;.    .lXMMMNx.   .xNO|", &
      "|:.   ;XMMMMMMd     :d.   ;KMMMMMMMMMMMX:    .dWMMWNl   .OO.      ;Xo   .kMMMMMW0o.    :XMNl.   ;0WMO|", &
      "|o'   .kMMMMMX:     :O,   .kMMMMMMMMMMNl      .kWMNXc   '0Wl     .kMd   .xMMMMMMMM0'   .kNl    cXMMMO|", &
      "|ko    ,0WMMXl.     ;Xd    'OWMMMNO0WNl   .'   .kWNK:   '0MK,    lNMx.   dMN0XWMMMK,   .kx.   :KWWWWk|", &
      "|OXl    .;c:'    '' ;XXl    .,cc:' .xl   .xK:   .OX0;   ,KMWk;,,c0MMx.   oWo..;cc:.    :Kc    .,,''''|", &
      "|OMNk,         .oKc ;XMNk,         ...  .dWMK,   ,kO;   ,KMMMWWWWMMMx.   o0,         .oXX;           |", &
      "|OMMMNOd:.    ;kNWl ;XMMMN0dlc::clxOkoookNMMMKdoodOXkoookNMMMMMMMMMMXxood0WKkoc:::cokXWMNkooooooooooc|", &
      "|OMMMMMMWx.    .,;. ;XMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMO|", &
      "|OMMMMMMMWKd;..   . ;XMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMO|", &
      "|OMMMMMMMMMMWX0kdd; ;XMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMO|", &
      "|OMMMMMMMMMMMMMMMWl ;XMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMWWWMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMO|", &
      "|OMMMMMMMMMMMMMMMWl ;XMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMKo;,;o0WMMMMMMMN00NMMMMMMMMMMMMMMMMMO|", &
      "|OMMMMMMMMMMMMMMMWl ;XMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMWk'.,:' .lXMMMMMNl..cXMMMMMMMMMMKodXMMO|", &
      "|OMMMMMMMMMMMMMMMWl ;XMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMWk..xNMNx' ,OWMMWd... cXMMMMMMMMMx..OMMO|", &
      "|OMMMMMMMMMMMMMMMWl ;XMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM0,.xWMMMMK: .kWMK,,K0,.dWMMMMMMMMx..OMMO|", &
      "|OMMMMMMMMMMMMMMMWl ;XMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMNc lNMMMMMMXl..xWx.oWMO.'0MMMMMMMMx..OMMO|", &
      "|OMMMMMMMMMMMMMMMWl ;XMMMMMMMMMMMMMMMMMMMMMMMMMMWNNNNNWMMMMMk.,0MMMMMMMMNo..d:'0MMWo.cNMMMMMMMx..OMMO|", &
      "|OMMMMMMMMMMMMMMMWl ;XMMMMMMMMMMMMMMMMMMMMMMNOo:;;;;;,,:okXX:.dWMMMMMMMMMNo...:NMMMK,.kMMMMMMMx..OMMO|", &
      "|OMMMMMMMMMMMMMMMWc ;XMMMMMMMMMMMMMMMMMMMWKo'.,ok0000kdc'..'.,KMMMMMMMMMMMNo  oWMMMMx.cNMMMMMMx..OMMO|", &
      "|OMMMMMMMMMMMMWKKNl ;XMMMMMMMMMMMMMMMMMW0c..cONMMMMMMMMMNO:. .cONMMMMMMMMMMX; 'OMMMMX;.kMMMMMMx..OMMO|", &
      "|OMMMMMMMMMMMM0,'kl ;XMMMMMMMMMMMMMMMW0c..c0WMMMMMMMMMMMMWo.';. ;kNMMMMMMMMK,  ,0MMMMx.cNMMMMMx..OMMO|", &
      "|OMMMMMMMMMMMMO..xl ;XMMMMMMMMMMMMMNk:..l0WMMMMMMMMMMMMMM0'.OWO:. ,kNMMMMMMx.:; ;KMMMX;.kMMMMMx..OMMO|", &
      "|OMMMMMMMMMMMMO..xl ;XMMMMMMMMWX0xc'.;dKWMMMMMMMMMMMMMMMNc.lNMMWO:. ;kNMMMWl.O0, :XMMMx.cNMMMMx..OMMO|", &
      "|OMMMMMMMMMMMMO..xl .lollllc:;;,,;cd0WMMMMMMMMMMMMMMMMMWd.;KMMMMMWO:. ;kNMK,;XMO' lNMMX;.kMMMMx..OMMO|", &
      "|OMMMMMMMMMMMMO..xl .cllloodxkOKXWMMMMMMMMMMWXNWMMMMMMWx..OMMMMMMMMNk, .:Oo.dWMMk..oNMWx.cNMMMx..OMMO|", &
      "|OMMMMMMMMMMMMO..xl ;XMMMMMMMMMMMMMMMMMMMMMMk,:KMMMMMXo..kWMMMMMMMMWx'.. ..,KMMMWx..dWMX;.kMMMx..OMMO|", &
      "|OMMMMMMMMMMMMO..xl ;XMMMMMMMMMMMMMMMMMMMMMMx..kK0xdc..;0WMMMMMMMMMWl ,x,  .lKMMMWx..xWWx.:NMMx..OMMO|", &
      "|OMMMMMMMMW0dOO'.xl ;XMMMMMMMMMMMMMMMMMMMMMMx. ......:kNMMMMMMMMMMMWl  . .;' .oXMMWd..xWX;.xMMx..OMMO|", &
      "|OMMMMMMMMWl :k'.xl ;XMMMMMMMMMMMMMMMMMMMMMMx..:oxOKXWMMMMMMMMMMMMMWl  .:OWXo. .oXMNo..dNx.;XMx..OMMO|", &
      "|OMMMMMMMMWl :k'.l; ;XMMMMMMMMMMMMMMMMMMMMMMx.'0MMMMMMMMMMMMMMMMMMMWl ,0WMMMWKo. .oKNd..oKc.oWx..OMMO|", &
      "|OMMMMMMMMWl :k'    ;XMMMMMMMMMMMMMMMMMMMMMMx.'0MMMMMMMMMMMMMMMMMMMWl ;XMMMMMMMXd' .cxd. ;l..xx..OMMO|", &
      "|OMMMMMMMMWl :k'    ;XMMMMMMMMMMMMMMMMMMMMMMx.'0MMMMMMMMMMMMMMMMMMMWl ;XMMMMMMMMMNk:..',. .. .' .OMMO|", &
      "|OMMMMMMMMWl :k'    ;XMMMMMMMMMMMMMMMMMMMMMMx.'0MMMMMMMMMMMMMMMMMMMWl ;XMMMMMMMMMMMWKxl,.       .OMMO|", &
      "|OMMMMMMMMWl :k'    'xXMMMMMMMMMMMMMMMMMMMMMx.'OMMMMMMMMMMMMMMMMMMMWl ;XMMMMMMMMMMMMMMMWX0kxdo, .OMMO|", &
      "|OMMMMMMMMWc :k'      oWMMMMMMMMMMMMMMMMMMMMx.'OMMMMMMMMMMMMMMMMMMMWl ;XMMMMMMMMMMMMMMMMMMMMMMx..OMMO|", &
      "|OMMMMMMMNO; :x.      lWMMMMMMMMMMMMMMMMMMMMx.'0MMMMMMMMMMMMMMMMMMMWl ;XMMMMMMMMMMMMMMMMMMMMMMx..xNWO|", &
      "|OMMMMMN0c.  ..       ,xXMMMMMNK0KNXxdKOcxX0c .OWWMMMMMMMMMMMMMMMMMWl ,ONMMMMMMMMMMMMMMMMMMMMMx. .cXO|", &
      "-----------------------------------------------------------------------------------------------------"
end subroutine qcxms2_logo

