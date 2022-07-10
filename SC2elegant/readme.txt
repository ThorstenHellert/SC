See also: Penn et al., "Comparisons Between AT and Elegant Tracking", (https://jacow.org/ipac2021/papers/MOPAB119.pdf)

The SC2elegant function takes an AT lattice stored in memory and
converts it to an elegant lattice file, along with separate files
for multipole fields as needed.

Elegant allows for multiple ways to represent a given magnet.
This code chooses a single format, and chooses elements with
symplectic methods such as CSBEND, KQUAD, and KSEXT.
CSBEND elements will have multiple fields expressed in terms of Fn and Gn
components.  The reference radius is hard wired to be 1 mm.
Apertures are included as needed, if they are present in the
AT definitions.

The code relies on metadata provided by the SC toolkit to
obtain magnet rolls and displacements.  The "malign" parameter
determines the specific type of error, and in the current form is fixed
to indicate body-centered misalignments.  When SC metadata is not
present, approximations are used for the misalignments and for large
errors they may not agree well.

The "group" for each element will remain the same after translation;
however, the name for copies of the same element is adjusted to have
a unique name for each magnet.  Thus, scripts in elegant that rely on
selection by element name will have to rely on wildcards.

Because the raw AT data structures are used, all modifications to a beamline
will be captured in this beamline.  However, there will be no procedural
definition of the lattice - identical copies of sections will each be
defined separately in the elegant file.

More runtime aspects of elegant, such as MALIGN elements
or "watch" elements will have to be inserted in the
lattice file either by hand or as part of the run script.  

One major aspect missing from the conversion process is the ability
to direct the nominal orbit out of the original plane, for example
to model vertical doglegs.  This capability will be added in the future.

The converter also does not handle BndStrMPoleSymplectic passmethods,
which should be converted into a CCBEND element in elegant,
Corrector, Matrix, Wiggler, or Thin passmethods.  Some of these
does have corresponding elements in elegant, so they could be added.

Comparisons between the original AT lattice and converted elegant
lattice can be perfomed directly.  However, because they codes have
different processes for dealing with closed orbits and phase offsets,
a few pre-adjustments to an AT ring will make the comparison process
simpler.

First, for a single RF cavity, it is better to adjust the TimeLag
value so that the calculated closed orbit has t=0 at the entrance to the
RF cavity.  This allows the TimeLag value to be directly converted into
an RF phase in elegant.

Second, it is helpful to adjust the RF frequency so that the closed
orbit in AT has 0 energy offset at some convenient location, either
the RF cavity or the s=0 position in the ring.

Example Matlab code for simple case with a single RF cavity:

circumf=findspos(RING,length(RING)+1)
rfLoc=findcells(RING,'PassMethod','RFCavityPass');
if length(rfLoc)>0
% adjust frequency to stay on momentum if wanted
   freq_temp=RING{rfLoc(1)}.Frequency;
   orb_temp=findorbit6(RING,[1]);
   orb_temp(5)=0; % track for nominal momentum instead
   c1turn=atpass(RING,orb_temp,1,2,[1]);
   RING{rfLoc(1)}.Frequency=freq_temp*circumf/(circumf+c1turn(6,2)-c1turn(6,1));
% adjust TimeLag to set closed orbit time to 0 at the RF cavity
   orb_temp=findorbit6(RING,[rfLoc(1)]);
   RING{rfLoc(1)}.TimeLag=RING{rfLoc(1)}.TimeLag-orb_temp(6);
end
