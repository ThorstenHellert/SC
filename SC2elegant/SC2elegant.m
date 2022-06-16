function SC2elegant( Lattice, LatticeName, FileName )
%
% Gregory Penn, 10 Oct 2021
%
% Makes a basic elegant input file out of an AT lattice data structure.
% Functions to handle specific elements are at the end of the file.
% Note:
%   Lattice is a Matlab object for the lattice (a list of structures)
%   LatticeName is a string which is what the lattice will be called
%   FileName is a string which is the name of the output file
%   for magnetic elements,
%     elegant ISR flag is being left at default value of 0
%   polyB is 'normal' components, polyA is 'skew' components
%   the 'tilt' parameter in elegant means rotate given fields by tilt
%      in the counterclockwise direction
%
   clear('elementlist')
   clear('newnamelist')
   clear('aperturelist')
   ndigits=3;
   xnom=0.001; % nominal radius for multipole fields (m)
   tiny=1.e-12;
   malign=2;   % rolls  0: old; 1: entrance-centered 2: body-centered
   refcor=1;   % 1 to use REFERENCE_CORRECTION; 0 not to use it
   energy=0;
% look for first value of energy
   for i=1:length(Lattice)
      if and(energy==0,isfield(Lattice{i},'Energy'))
         energy=Lattice{i}.Energy;
      end
   end
   msg=sprintf('Energy for lattice is %.15g\n',energy);
   disp(msg);
% dictionary of original names:  
%   how many times appears, how many times used to generate a new name, 
%   and how many digits used in element name
   oldnamedict=containers.Map();
   newnamelist=[];
   for i=1:length(Lattice)
      elementlist{i}=ConstrainName(Lattice{i}.FamName);
      oldnamedict(elementlist{i})=[0,0,ndigits];
   end
% count how many iterations of each element
   for i=1:length(Lattice)
      counts=oldnamedict(elementlist{i});
      oldnamedict(elementlist{i})=[counts(1)+1,counts(2),counts(3)];
   end
% define new namelist
   for i=1:length(Lattice)
      oldname=elementlist{i};
      counts=oldnamedict(oldname);
      ndigcur=counts(3);
      formatstr=['%0' num2str(ndigcur) '.f'];
      if counts(1)==1
         newnamelist{i}=oldname;
      else
         counts(2)=counts(2)+1;
         oldnamedict(oldname)=counts;
         newname=[oldname '_' num2str(counts(2),formatstr)];
         while (ismember(newname,newnamelist))
            counts(2)=counts(2)+1;
            if counts(2)>=10^ndigcur
	       ndigcur=ndigcur+1;
               counts(3)=ndigcur;
               formatstr=['%0' num2str(ndigcur) '.f'];
            end
            newname=[oldname '_' num2str(counts(2),formatstr)];
         end
         oldnamedict(oldname)=counts;
         newnamelist{i}=newname;
      end
   end
% prep for apertures as needed
   aperturelist=cell(1,length(Lattice));
% requires actively setting needs_aperture to 1 to place MAXAMP element
   needs_aperture=zeros(1,length(Lattice));
   % for when no aperture is defined
   aperture.Type='N';
   aperture.Values=[];
   aperture.Name='';
   aperture_active=aperture;
   for i=1:length(Lattice)
      if isfield(Lattice{i},'RApertures')
         aperture.Type='R';
         aperture.Values=Lattice{i}.RApertures;
         needs_aperture(i)=1;
      elseif isfield(Lattice{i},'EApertures')
         aperture.Type='E';
         aperture.Values=Lattice{i}.EApertures;
         needs_aperture(i)=1;
      end
      % skip if there is no actual change in the aperture
      if and(isequal(aperture.Type,aperture_active.Type), ...
             isequal(aperture.Values,aperture_active.Values))
         needs_aperture(i)=0;
      end
      if needs_aperture(i)==1
         aperture.Name=['APER_' newnamelist{i}];
         count=0;
         while (ismember(aperture.Name,newnamelist))
            count=count+1;
            aperture.Name=['APER_' newnamelist{i} '_' num2str(count)];
         end
         aperturelist{i}=aperture;
         aperture_active=aperture;
      end
   end
%
% write elegant input file
   fid=fopen(FileName,'wt+');
   fprintf(fid, '! elegant input file based on AT beamline\n');
   fprintf(fid,'! energy\t%.15g\n',energy);
   fprintf(fid, '!\n');
   fprintf(fid, '! defining the standard individual elements\n');
% define each standard element and write to the file
   for i=1:length(Lattice)
      switch (Lattice{i}.PassMethod)
      case 'IdentityPass'
         MakeMarker(Lattice{i},newnamelist{i},fid);
      case 'DriftPass'
         MakeDrift(Lattice{i},newnamelist{i},fid);
      case 'CavityPass'
         MakeRF(Lattice{i},newnamelist{i},fid);
      case 'RFCavityPass'
         MakeRF(Lattice{i},newnamelist{i},fid);
      case 'StrMPoleSymplectic4Pass'
         radloss=0;
         switch Lattice{i}.Class;
         case 'Quadrupole'
            MakeQuad(Lattice{i},newnamelist{i},elementlist{i},fid,radloss);
         case 'Sextupole'
            MakeSext(Lattice{i},newnamelist{i},elementlist{i},fid,radloss);
         case 'Octupole'
            MakeOct(Lattice{i},newnamelist{i},elementlist{i},fid,radloss);
         otherwise
            msg=sprintf('Warning:  unknown type %s in element %s; treating as quadrupole\n',Lattice{i}.Class,newnamelist{i});
            disp(msg);
            MakeQuad(Lattice{i},newnamelist{i},elementlist{i},fid,radloss);
         end 
      case 'StrMPoleSymplectic4RadPass'
         radloss=1;
         switch Lattice{i}.Class;
         case 'Quadrupole'
            MakeQuad(Lattice{i},newnamelist{i},elementlist{i},fid,radloss);
         case 'Sextupole'
            MakeSext(Lattice{i},newnamelist{i},elementlist{i},fid,radloss);
         case 'Octupole'
            MakeOct(Lattice{i},newnamelist{i},elementlist{i},fid,radloss);
         otherwise
            msg=sprintf('Warning:  unknown type %s in element %s; treating as quadrupole\n',Lattice{i}.Class,newnamelist{i});
            disp(msg);
            MakeQuad(Lattice{i},newnamelist{i},elementlist{i},fid,radloss);
         end
      case 'BndMPoleSymplectic4Pass'
         radloss=0;
         MakeDipole(Lattice{i},newnamelist{i},elementlist{i},fid,radloss);
      case 'BndMPoleSymplectic4RadPass'
         radloss=1;
         MakeDipole(Lattice{i},newnamelist{i},elementlist{i},fid,radloss);
      otherwise
         msg=sprintf('Warning:  unknown method %s for element %s; treating as drift\n',Lattice{i}.PassMethod,newnamelist{i});
         disp(msg);
         MakeDrift(Lattice{i},newnamelist{i},fid);
      end
   end
%
% write each aperture element as needed to the file
   fprintf(fid, '!\n');
   fprintf(fid, '! defining the apertures as needed\n');
   for i=1:length(Lattice)
      if and(needs_aperture(i)==1,length(aperturelist{i})>0)
         MakeAperture(aperturelist{i},fid);
      end
   end
%
% make the ring element
   fprintf(fid, '!\n');
   fprintf(fid, '! defining the full beamline\n');
   fprintf(fid, '%s: LINE=(&\n',LatticeName);
   for i=1:length(Lattice)
      endstr=',&';
      if i==length(Lattice)
         endstr=')';
      end  
      if needs_aperture(i)==0
         fprintf(fid,'%s%s\n',newnamelist{i},endstr);
      else
         fprintf(fid,'%s,%s%s\n',aperturelist{i}.Name,newnamelist{i},endstr);
      end
   end
   fclose(fid);
%
% helper functions
%
% Marker element
   function MakeMarker(t_element,t_name,t_id);
      fprintf(t_id,'%s: MARK\n',t_name);
   end
%
% Drift element
   function MakeDrift(t_element,t_name,t_id);
      is_col=0;
      % sometimes when there are apertures, may be treated as a collimator 
      if isfield(t_element,'RApertures')
         ApertureVals=t_element.RApertures;
         xmin=ApertureVals(1);
         xmax=ApertureVals(2);
         ymin=ApertureVals(3);
         ymax=ApertureVals(4);
         if (xmax==-xmin)
            xm=xmax;
            dx=0;
         else
	    xm=0.5*(xmax-xmin);
            dx=0.5*(xmax+xmin);
         end
         if (ymax==-ymin)
            ym=ymax;
            dy=0;
         else
	    ym=0.5*(ymax-ymin);
            dy=0.5*(ymax+ymin);
         end
         if t_element.Length==0
            is_col=1;
         elseif dx~=0
            is_col=1;
         elseif dy~=0
            is_col=1;
         end
         if and(dx==0,dy==0)
            colstring=sprintf('%s: RCOL,L=%.15g,X_MAX=%.15g,Y_MAX=%.15g',t_name,t_element.Length,xm,ym);
         else
            colstring=sprintf('%s: RCOL,L=%.15g,X_MAX=%.15g,Y_MAX=%.15g,DX=%.15g,DY=%.15g',t_name,t_element.Length,xm,ym,dx,dy);
         end
      elseif isfield(t_element,'EApertures')
         ApertureVals=t_element.EApertures;
         xmax=ApertureVals(1);
         ymax=ApertureVals(2);
         if t_element.Length==0
            is_col=1;
            colstring=sprintf('%s: ECOL,L=%.15g,X_MAX=%.15g,Y_MAX=%.15g',t_name,t_element.Length,xmax,ymax);
         end
      end
      if is_col==0
         fprintf(t_id,'%s: EDRIFT,L=%.15g\n',t_name,t_element.Length);
      else
         fprintf(t_id,'%s\n',colstring);
      end
   end
%
% RF cavity element - is TimeLag how phase is set?
   function MakeRF(t_element,t_name,t_id);
      % phase is in degrees
      phase=180.-360.*t_element.Frequency*t_element.TimeLag/299792458.;
      % multiply frequency by v/c, AT assumes v=c
      bg0=energy/(0.51099906*1.e6);
      b0=bg0/sqrt(1+bg0^2);
      adjFreq=t_element.Frequency*b0;
%      fprintf(t_id,'%s: RFCA,L=%.15g,FREQ=%.20g,VOLT=%.15g,PHASE=%.20g,CHANGE_T=0\n',t_name,t_element.Length,t_element.Frequency,t_element.Voltage,phase);
      fprintf(t_id,'%s: RFCA,L=%.15g,FREQ=%.20g,VOLT=%.15g,PHASE=%.20g,CHANGE_T=0\n',t_name,t_element.Length,adjFreq,t_element.Voltage,phase);
   end
%
% Dipole
   function MakeDipole(t_element,t_name,t_group,t_id,t_radloss);
      angle=t_element.BendingAngle;
      e1=t_element.EntranceAngle;
      e2=t_element.ExitAngle;
      % inverse of nominal bending radius (so nominal version of k0)
      % same sign as bending angle
      invrhoN=angle/t_element.Length;
      polyB=t_element.PolynomB;
      polyA=t_element.PolynomA;
      % note that polyB is the normal component, polyA is skew
      % separate out vertical kick
%      dipole_vkick=0;
%      dipole_hkick=0;
      dipole_vkick=polyA(1)*t_element.Length;
      polyA(1)=0;
      dipole_hkick=-polyB(1)*t_element.Length;
      polyB(1)=0;
      n_fields=length(polyB);
      if n_fields>9
         n_fields=9;
      end
      % assuming tilt=0 for now;  will obtain from rotation matrices
      tilt=0;
      nom_k0=invrhoN;
      nom_k1=0;
      % seems that CANNOT use K1 field at same time as using F1
      % so there is no point in trying to specify nominal value
      %   should check with Michael Borland
%      if isfield(t_element,'K')
%         nom_k1=t_element.K;
%      end
      nom_polyB=[0];
      nom_polyA=[0];
      etilt=0;
      % read nominal fields if available
      %  though not using nominal dipole components
      %    does not help to distinguish etilt from vertical corrector
      %    also intentional tilt in AT must be done using rotation elements
      if isfield(t_element,'NomPolynomB')
         nom_polyB=t_element.NomPolynomB;
         nom_polyA=t_element.NomPolynomA;
      end
      kickB=polyB(1);
      kickA=polyA(1);
      % difference in full fields modelled as error (etilt, fse_dipole)
      etiltA=-atan2(kickA,kickB+nom_k0);
      k0=sqrt((kickB+nom_k0)^2+kickA^2);
      if (abs(etiltA)>pi/2)
         etiltA=-atan2(-kickA,-(kickB+nom_k0));
         k0=-k0;
      end
      etilt=etilt+etiltA;
      str_offset=HandleOffsets(t_element,t_name);
      % currently not handling intentional vertical component of reference
      [str_roll, errtilt]=HandleRolls(t_element,t_name,1);
      fse_dipole=k0/nom_k0-1;
      % transform higher order fields by rotating counterclockwise by etilt
      % elegant will rotate them in the opposite direction      
      if n_fields>=2
         fields_norm=zeros(n_fields-1,1);
         fields_skew=zeros(n_fields-1,1);
         for j=2:n_fields
            fields_norm(j-1)=polyB(j)*cos(j*etilt)-polyA(j)*sin(j*etilt);
            fields_skew(j-1)=-(polyA(j)*cos(j*etilt)+polyB(j)*sin(j*etilt));
         end
         % subtract nominal k1
         fields_norm(1)=fields_norm(1)-nom_k1;
         % rescale for elegant
         for j=2:n_fields
            % definitely nom_k0 works better than k0 for quad and sextupole
            % fields when using FSE_DIPOLE, but for now using FSE so ok
            if (j<2)
               fields_norm(j-1)=fields_norm(j-1)*xnom^(j-1)/nom_k0;
               fields_skew(j-1)=fields_skew(j-1)*xnom^(j-1)/nom_k0;
            else
               fields_norm(j-1)=fields_norm(j-1)*xnom^(j-1)/k0;
               fields_skew(j-1)=fields_skew(j-1)*xnom^(j-1)/k0;
            end
         end
      end
      str_errors='';
      if strcmp(str_roll,'')  % assume no 2nd etilt for now
	 if (fse_dipole~=0)
            t_str=sprintf('ETILT=%.15g,FSE_DIPOLE=%.15g,',etilt,fse_dipole);
            str_errors=[str_errors t_str];
         end
      elseif or((etilt~=0),(fse_dipole~=0))
         t_str=sprintf('FSE_DIPOLE=%.15g,',fse_dipole);
         str_errors=[str_errors t_str];
      end
      if n_fields>=2
         t_str=sprintf('XREFERENCE=%.15g,',xnom);
         str_errors=[str_errors t_str];
         for j=1:(n_fields-1)
            if fields_norm(j)~=0
              t_str=sprintf('F%d=%.15g,',j,fields_norm(j));
              str_errors=[str_errors t_str];
            end
            if fields_skew(j)~=0
              t_str=sprintf('G%d=%.15g,',j,fields_skew(j));
              str_errors=[str_errors t_str];
            end
         end
      end
% reference_correction used
      if and(dipole_vkick==0,dipole_hkick==0)
         fprintf(t_id,'%s: CSBEND,GROUP=%s,L=%.15g,ANGLE=%.15g,K1=%.15g,E1=%.15g,E2=%.15g,TILT=%.15g,%s%s%sEDGE1_EFFECTS=3,EDGE2_EFFECTS=3,REFERENCE_CORRECTION=%d,MALIGN_METHOD=%d,SYNCH_RAD=%d,N_SLICES=%d\n',t_name,t_group,t_element.Length,angle,nom_k1,e1,e2,tilt,str_offset,str_roll,str_errors,refcor,malign,t_radloss,t_element.NumIntSteps);
      else
         fprintf(t_id,'%s: CSBEND,GROUP=%s,L=%.15g,ANGLE=%.15g,K1=%.15g,E1=%.15g,E2=%.15g,XKICK=%.15g,YKICK=%.15g,TRACKING_MATRIX=3,TILT=%.15g,%s%s%sEDGE1_EFFECTS=3,EDGE2_EFFECTS=3,REFERENCE_CORRECTION=%d,MALIGN_METHOD=%d,SYNCH_RAD=%d,N_SLICES=%d\n',t_name,t_group,t_element.Length,angle,nom_k1,e1,e2,dipole_hkick,dipole_vkick,tilt,str_offset,str_roll,str_errors,refcor,malign,t_radloss,t_element.NumIntSteps);
      end
   end
%
% Quadrupole
   function MakeQuad(t_element,t_name,t_group,t_id,t_radloss);
      polyB=t_element.PolynomB;
      polyA=t_element.PolynomA;
      % note that polyB is the normal component, polyA is skew
      n_fields=length(polyB);
      nom_order=1;
      nom_i=nom_order+1;
      tilt=0;
      str_offset=HandleOffsets(t_element,t_name);
      % for elegant, will add roll to intentional tilt (e.g., skew)
      [str_roll, tilt]=HandleRolls(t_element,t_name,0);
      % check whether errors require auxiliary multipole file
      is_simple=1;
      if (n_fields>nom_i)
	if or(any(polyB((nom_i+1):n_fields)),any(polyA((nom_i+1):n_fields)))
           is_simple=0;
        end
      end
      if (is_simple)
         k1=sqrt(polyA(nom_i)^2+polyB(nom_i)^2);
         tiltA=-atan2(polyA(nom_i),polyB(nom_i))/nom_i;
         if (abs(tiltA)>pi/(2*nom_i))
            tiltA=-atan2(-polyA(nom_i),-polyB(nom_i))/nom_i;
            k1=-k1;
         end
         tilt=tilt+tiltA;
	 % convert k to elegant scaling (not relevant for k1)
         k1=k1*factorial(nom_order);
         % include dipole kicks if any
         str_kick=HandleKicks(t_element,tiltA);
         if strcmp(str_roll,'')  % there could be a tilt from the fields
                                 % assume there will not be 2 sources of tilt
	    str_roll=sprintf('TILT=%.15g,',tilt);
         else % need to specify geometric error model
            str_roll=[str_roll sprintf('MALIGN_METHOD=%d,',malign)];
         end
         fprintf(t_id,'%s: KQUAD,GROUP=%s,L=%.15g,K1=%.15g,%s%s%sSYNCH_RAD=%d,N_SLICES=%d\n',t_name,t_group,t_element.Length,k1,str_kick,str_offset,str_roll,t_radloss,t_element.NumIntSteps);
      else
         [k1 tiltA systematic]=CalcMultipole(t_element,nom_i,xnom);
         tilt=tilt+tiltA;
         % include dipole kicks if any
         str_kick=HandleKicks(t_element,tiltA);
         str_multi=[t_name '_multi.sdds'];
         if strcmp(str_roll,'')  % there could be a tilt from the fields
                                 % assume there will not be 2 sources of tilt
	    str_roll=sprintf('TILT=%.15g,',tilt);
         else % need to specify geometric error model
            str_roll=[str_roll sprintf('MALIGN_METHOD=%d,',malign)];
         end
         fprintf(t_id,'%s: KQUAD,GROUP=%s,L=%.15g,K1=%.15g,%s%s%sSYNCH_RAD=%d,N_SLICES=%d,SYSTEMATIC_MULTIPOLES="%s"\n',t_name,t_group,t_element.Length,k1,str_kick,str_offset,str_roll,t_radloss,t_element.NumIntSteps,str_multi);
         % write the systematics to file
         WriteMultipole(systematic,str_multi,xnom);
      end
   end
%
% Sextupole
   function MakeSext(t_element,t_name,t_group,t_id,t_radloss);
      polyB=t_element.PolynomB;
      polyA=t_element.PolynomA;
      % note that polyB is the normal component, polyA is skew
      n_fields=length(polyB);
      nom_order=2;
      nom_i=nom_order+1;
      tilt=0;
      str_offset=HandleOffsets(t_element,t_name);
      % for elegant, will add roll to intentional tilt (e.g., skew)
      [str_roll, tilt]=HandleRolls(t_element,t_name,0);
      % check whether errors require auxiliary multipole file
      %  ksext element uniquely can handle weak quad components 
      %  without external file
      is_simple=1;
      if (n_fields>nom_i)
	if or(any(polyB((nom_i+1):n_fields)),any(polyA((nom_i+1):n_fields)))
           is_simple=0;
        end
      end
      if is_simple
         k2=sqrt(polyA(nom_i)^2+polyB(nom_i)^2);
         tiltA=-atan2(polyA(nom_i),polyB(nom_i))/nom_i;
         if (abs(tiltA)>pi/(2*nom_i))
            tiltA=-atan2(-polyA(nom_i),-polyB(nom_i))/nom_i;
            k2=-k2;
         end
         tilt=tilt+tiltA;
         % convert k to elegant scaling
         k2=k2*factorial(nom_order);
         % include dipole kicks if any
         str_kick=HandleKicks(t_element,tiltA);
         % transform quadrupole fields by rotating counterclockwise by tilt
         k1=polyB(2)*cos(2*tiltA)-polyA(2)*sin(2*tiltA);
         j1=-(polyA(2)*cos(2*tiltA)+polyB(2)*sin(2*tiltA));
         if strcmp(str_roll,'')  % there could be a tilt from the fields
                                 % assume there will not be 2 sources of tilt
	    str_roll=sprintf('TILT=%.15g,',tilt);
         else % need to specify geometric error model
            str_roll=[str_roll sprintf('MALIGN_METHOD=%d,',malign)];
         end
         fprintf(t_id,'%s: KSEXT,GROUP=%s,L=%.15g,K2=%.15g,K1=%.15g,J1=%.15g,%s%s%sSYNCH_RAD=%d,N_SLICES=%d\n',t_name,t_group,t_element.Length,k2,k1,j1,str_kick,str_offset,str_roll,t_radloss,t_element.NumIntSteps);
      else
         [k2 tiltA systematic]=CalcMultipole(t_element,nom_i,xnom);
         % include dipole kicks if any
         str_kick=HandleKicks(t_element,tiltA);
         str_multi=[t_name '_multi.sdds'];
         if strcmp(str_roll,'')  % there could be a tilt from the fields
                                 % assume there will not be 2 sources of tilt
	    str_roll=sprintf('TILT=%.15g,',tilt);
         else % need to specify geometric error model
            str_roll=[str_roll sprintf('MALIGN_METHOD=%d,',malign)];
         end
         fprintf(t_id,'%s: KSEXT,GROUP=%s,L=%.15g,K2=%.15g,%s%s%sSYNCH_RAD=%d,N_SLICES=%d,SYSTEMATIC_MULTIPOLES="%s"\n',t_name,t_group,t_element.Length,k2,str_kick,str_offset,str_roll,t_radloss,t_element.NumIntSteps,str_multi);
         % write the systematics to file
         WriteMultipole(systematic,str_multi,xnom);
      end
   end
%
% Octupole
   function MakeOct(t_element,t_name,t_group,t_id,t_radloss);
      polyB=t_element.PolynomB;
      polyA=t_element.PolynomA;
      % note that polyB is the normal component, polyA is skew
      n_fields=length(polyB);
      nom_order=3;
      nom_i=nom_order+1;
      tilt=0;
      str_offset=HandleOffsets(t_element,t_name);
      % for elegant, will add roll to intentional tilt (e.g., skew)
      [str_roll, tilt]=HandleRolls(t_element,t_name,0);
      % dipole fields?
      if or(polyB(1)~=0,polyA(1)~=0)
         msg=sprintf('Warning:  octupole cannot have correction kicks, dipole fields will be ignored in element %s',t_name);
         disp(msg);
      end
      % check whether errors require auxiliary multipole file
      is_simple=1;
      if or(any(polyB(1:(nom_i-1))),any(polyA(1:(nom_i-1))))
         is_simple=0;
      end
      if (n_fields>nom_i)
	if or(any(polyB((nom_i+1):n_fields)),any(polyA((nom_i+1):n_fields)))
           is_simple=0;
        end
      end
      if is_simple
         k3=sqrt(polyA(nom_i)^2+polyB(nom_i)^2);
         tiltA=-atan2(polyA(nom_i),polyB(nom_i))/nom_i;
         if (abs(tilt)>pi/(2*nom_i))
            tiltA=-atan2(-polyA(nom_i),-polyB(nom_i))/nom_i;
            k1=-k1;
         end
         tilt=tilt+tiltA;
         % convert k to elegant scaling
         k3=k3*factorial(nom_order);
         if strcmp(str_roll,'')  % there could be a tilt from the fields
                                 % assume there will not be 2 sources of tilt
	    str_roll=sprintf('TILT=%.15g,',tilt);
         else % need to specify geometric error model
            str_roll=[str_roll sprintf('MALIGN_METHOD=%d,',malign)];
         end
         fprintf(t_id,'%s: KOCT,GROUP=%s,L=%.15g,K3=%.15g,%s%sSYNCH_RAD=%d,N_SLICES=%d\n',t_name,t_group,t_element.Length,k3,str_offset,str_roll,t_radloss,t_element.NumIntSteps);
      else
         [k3 tiltA systematic]=CalcMultipole(t_element,nom_i,xnom);
         str_multi=[t_name '_multi.sdds'];
         if strcmp(str_roll,'')  % there could be a tilt from the fields
                                 % assume there will not be 2 sources of tilt
	    str_roll=sprintf('TILT=%.15g,',tilt);
         else % need to specify geometric error model
            str_roll=[str_roll sprintf('MALIGN_METHOD=%d,',malign)];
         end
         fprintf(t_id,'%s: KOCT,GROUP=%s,L=%.15g,K3=%.15g,%s%sSYNCH_RAD=%d,N_SLICES=%d,SYSTEMATIC_MULTIPOLES="%s"\n',t_name,t_group,t_element.Length,k3,str_offset,str_roll,t_radloss,t_element.NumIntSteps,str_multi);
         % write the systematics to file
         WriteMultipole(systematic,str_multi,xnom);
      end
   end
%
% Aperture element - has to be separated from the originating AT element
   function MakeAperture(t_aper,t_id);
      switch t_aper.Type
      case 'R'
         xmin=t_aper.Values(1);
         xmax=t_aper.Values(2);
         ymin=t_aper.Values(3);
         ymax=t_aper.Values(4);
         if and(xmax==-xmin,ymax==-ymin)
            fprintf(t_id,'%s: MAXAMP,X_MAX=%.15g,Y_MAX=%.15g,ELLIPTICAL=0\n',t_aper.Name,xmax,ymax);
         else
            openside='';
            xm=xmax;
            ym=ymax;
            if (xmax>-xmin)
               openside=[openside '+x'];
               xm=-xmin;
            elseif (xmax<-xmin)
               openside=[openside '-x'];
               xm=xmax;
            end
            if (ymax>-ymin)
               openside=[openside '+y'];
               ym=-ymin;
            elseif (ymax<-ymin)
               openside=[openside '-y'];
               ym=ymax;
            end
            fprintf(t_id,'%s: MAXAMP,X_MAX=%.15g,Y_MAX=%.15g,ELLIPTICAL=0,OPEN_SIDE=%s\n',t_aper.Name,xm,ym,openside);
         end
      case 'E'
         fprintf(t_id,'%s: MAXAMP,X_MAX=%.15g,Y_MAX=%.15g,ELLIPTICAL=1\n',t_aper.Name,t_aper.Values(1),t_aper.Values(2));
      otherwise
         msg=sprintf('Warning:  unknown aperture type for %s; treating as MARK instead\n',t_aper.Name);
         disp(msg);
         fprintf(t_id,'%s: MARK\n',t_aper.Name);
      end
   end
%
% Treat offsets
   function str_out = HandleOffsets(t_element,t_name);
      dx=0;
      dy=0;
      dz=0;
      if or(isfield(t_element,'MagnetOffset'),isfield(t_element,'SupportOffset'))
         if isfield(t_element,'MagnetOffset')
            magoff=threeVec(t_element.MagnetOffset);
            dx=dx+magoff(1);
            dy=dy+magoff(2);
            dz=dz+magoff(3);
         end
         if isfield(t_element,'SupportOffset')
            supoff=threeVec(t_element.SupportOffset);
            dx=dx+supoff(1);
            dy=dy+supoff(2);
            dz=dz+supoff(3);
         end
      else
      % rely on T1,T2 if no other fields giving information on displacements
         if isfield(t_element,'T1')
	    match=true;
	    dx=t_element.T1(1);
	    dy=t_element.T1(3);
	    dz=t_element.T1(5);
            if max(abs(t_element.T1+t_element.T2))>tiny
	       match=false;
	       dx=0.5*(t_element.T1(1)-t_element.T2(1));
	       dy=0.5*(t_element.T1(3)-t_element.T2(3));
	       dz=0.5*(t_element.T1(5)-t_element.T2(5));
               msg=sprintf('Warning: incoming and outgoing displacements do not match for element %s; using the average value\n',t_name);
               disp(msg);
	    end               
         end
      end
      if any([dx dy dz])
         str_out=sprintf('DX=%.15g,DY=%.15g,DZ=%.15g,',dx,dy,dz);
      else
         str_out='';
      end
   end
%
% Treat rolls
   function [str_out, errtilt] = HandleRolls(t_element,t_name,is_bend);
      ax=0;
      ay=0;
      az=0;
      if or(isfield(t_element,'MagnetRoll'),isfield(t_element,'SupportRoll'))
         if isfield(t_element,'MagnetRoll')
            magroll=threeVec(t_element.MagnetRoll);
	    ax=ax+magroll(2);
            ay=ay+magroll(3);
            az=az+magroll(1);
         end
         if isfield(t_element,'SupportRoll')
            suproll=threeVec(t_element.SupportRoll);
            ax=ax+suproll(2);
            ay=ay+suproll(3);
            az=az+suproll(1);
         end
      else
      % rely on R1,R2 if no other fields giving information on displacements
      % should do this better
         if isfield(t_element,'R1')
            msg=sprintf('Warning: trying to match rotation matrix without information about MagnetRoll or SupportRoll in element %s\n',t_name);
            disp(msg);
	    ax=-atan(t_element.R1(2,3)/t_element.R1(3,3));
	    ay=-asin(t_element.R1(1,3));
            az=-atan(t_element.R1(1,2)/t_element.R1(1,1));
         end
      end
      errtilt=az;
      str_out='';
      if (az~=0)
         if is_bend
            str_out=[str_out sprintf('ETILT=%.15g,',az)];
         else
            str_out=[str_out sprintf('TILT=%.15g,',az)];
         end
      end
      if (ax~=0)
         if is_bend
            str_out=[str_out sprintf('EPITCH=%.15g,',ax)];
         else
            str_out=[str_out sprintf('PITCH=%.15g,',ax)];
         end
      end
      if (ay~=0)
         if is_bend
            str_out=[str_out sprintf('EYAW=%.15g,',ay)];
         else
            str_out=[str_out sprintf('YAW=%.15g,',ay)];
         end
      end
   end
%
% Treat dipole kicks in straight elements
   function str_out = HandleKicks(t_element,t_tilt);
      polyB=t_element.PolynomB;
      polyA=t_element.PolynomA;
      % note that polyB is the normal component, polyA is skew
      if or(polyB(1)~=0,polyA(1)~=0)
         kick_norm=polyB(1)*cos(t_tilt)-polyA(1)*sin(t_tilt);
         kick_skew=-(polyA(1)*cos(t_tilt)+polyB(1)*sin(t_tilt));
         hkick=-kick_norm*t_element.Length;
         vkick=-kick_skew*t_element.Length;
         str_out=sprintf('HKICK=%.15g,VKICK=%.15g,',hkick,vkick);
      else
         str_out='';
      end
   end
%
% Calculate "systematic" multipoles (excludes dipole field)
   function [t_k t_tilt t_systematic] = CalcMultipole(t_element,t_i,t_x);
      polyB=t_element.PolynomB;
      polyA=t_element.PolynomA;
      n_fields=length(polyB);
      t_k=sqrt(polyB(t_i)^2+polyA(t_i)^2);
      % start with nominal fields if available
      if (t_i==1)
         t_tilt=0;
      elseif isfield(t_element,'NomPolynomB')
         nom_norm=t_element.NomPolynomB(t_i);
         nom_skew=t_element.NomPolynomA(t_i);
         t_k=sqrt(nom_norm^2+nom_skew^2);
         t_tilt=-atan2(nom_skew,nom_norm)/t_i;
         if (abs(t_tilt)>pi/(2*t_i))
            t_tilt=-atan2(-nom_skew,-nom_norm)/t_i;
            t_k=-t_k;
         end
         % subtract the nominal values for the primary field
         polyB(t_i)=polyB(t_i)-nom_norm;
         polyA(t_i)=polyA(t_i)-nom_skew;
      else
         % rotate entire field of nominal index t_i
	 t_k=sqrt(polyB(t_i)^2+polyA(t_i)^2);
         t_tilt=-atan2(polyA(t_i),polyB(t_i))/t_i;
         if (abs(t_tilt)>pi/(2*t_i))
            t_tilt=-atan2(-polyA(t_i),-polyB(t_i))/t_i;
            t_k=-t_k;
         end
         % remove the primary field, treat other orders as systematics
         polyB(t_i)=0;
         polyA(t_i)=0;
      end
      % scaling: multiply field of order t_i-1 by factorial(t_i-1)
      t_k=t_k*factorial(t_i-1);
      % rotate remaining fields counterclockwise by t_tilt
      %  elegant will rotate them in opposite direction
      %  also convert to elegant format
      %    looks like polyA has opposite sign to elegant skew
      t_systematic=zeros(n_fields-1,2);
      for j=2:n_fields
         t_systematic(j-1,1)=polyB(j)*cos(j*t_tilt)-polyA(j)*sin(j*t_tilt);
         t_systematic(j-1,2)=-(polyA(j)*cos(j*t_tilt)+polyB(j)*sin(j*t_tilt));
      end  
      for j=2:n_fields
	 t_systematic(j-1,1)=t_systematic(j-1,1)*(t_x)^(j-t_i)/t_k*factorial(t_i-1);
         t_systematic(j-1,2)=t_systematic(j-1,2)*(t_x)^(j-t_i)/t_k*factorial(t_i-1);
      end
   end
%
% Write "systematic" multipole info to file
   function WriteMultipole(multipoles,t_name,t_x);
      fileid=fopen(t_name,'wt');
      [n_fields ncol]=size(multipoles);
      fprintf(fileid,'SDDS1\n');
%      fprintf(fileid,'&parameter name=referenceRadius, units=m, type=double, fixed_value=%g &end\n',t_x);
      fprintf(fileid,'&parameter name=referenceRadius, units=m, type=double &end\n',t_x);
      fprintf(fileid,'&column name=order, type=long &end\n');
      fprintf(fileid,'&column name=normal, type=double &end\n');
      fprintf(fileid,'&column name=skew, type=double &end\n');
      fprintf(fileid,'&data mode=ascii &end\n');
      fprintf(fileid,'%.15g\n',t_x);
      fprintf(fileid,'\t%d\n',n_fields);
      for j=1:n_fields
         fprintf(fileid,'%d\t%.15g\t%.15g\n',j,multipoles(j,1),multipoles(j,2));
      end
%      fprintf(fileid,'&data mode=binary, column_major_order=1 &end\n');
%      fwrite(fileid, t_x, 'double');
%      fwrite(fileid, n_fields, 'int32');
%      fwrite(fileid, (1:n_fields), 'long');
%      fwrite(fileid, real(transpose(multipoles(:,1))), 'double');
%      fwrite(fileid, real(transpose(multipoles(:,2))), 'double');
      fclose(fileid);      
   end
%
% Ensure that element names are consistent with elegant
   function str_name = ConstrainName(t_name);
      str_name=t_name;
      if ~isletter(t_name(1))
         str_name=['A' str_name];
      end
      % semicolon is not actually forbidden but is a pain
      %    requires quotes around every name that contains it
      %    assign a different replacement for each unwanted character
      for j=1:length(str_name)
         switch str_name(j)
         case ':'
            str_name(j)='|';
         case '#'
            str_name(j)='$';
         case '*'
            str_name(j)='&';
         case '!'
            str_name(j)='?';
         case char(34)
            str_name(j)='^';
         case char(39)
            str_name(j)='@';
         case char(44)
            str_name(j)='.';
         case '('
            str_name(j)='[';
         case ')'
            str_name(j)=']';
         end
      end
   end
%
% function to ensure x,y,z coordinates are produced
% even when fewer or more are given; assumes x,y,z in that order
   function outvec = threeVec(t_vec);
      outvec=zeros(1,3);
      for i=1:min(3,length(t_vec))
         outvec(i)=t_vec(i);
      end
   end

end
