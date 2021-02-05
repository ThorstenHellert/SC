function [T1,T2,R1,R2] = SCgetTransformation(dx,dy,dz,ax,ay,az,magTheta,magLength,varargin)
% SCgetTransformation
% ===================
%
% NAME
% ----
% SCgetTransformation - calculates the translation and rotation elements for given misalignments
%
% SYNOPSIS
% --------
% `[T1,T2,R1,R2] = SCgetTransformation(dx,dy,dz,ax,ay,az,magTheta,magLength [,options])`
%
%
% DESCRIPTION
% -----------
% This function calculates AT's linear rotation and translation elements `T1`, `T2`, and
% `R1`, `R2` for magnets with aribtrary horizontal, vertical and longitudinal offsets `[dx,dy,dz]`
% and roll `az` (roll around z-axis), pitch `ax` (roll around x-axis) and yaw `ay` (roll around
% y-axis), see Tech-Note XXX for more details. By default all errors are defined in the midpoint of
% a straight line connecting the ideal entrance and exit point of the magnet.
%
% INPUT
% -----
% `dx`::
%	Horizontal offset.
% `dy`::
%	Vertical offset.
% `dz`::
%	Longitudinal offset.
% `ax`::
%	Roll around x-axis ('pitch')
% `ay`::
%	Roll around y-axis ('yaw')
% `az`::
%	Roll around z-axis ('roll')
% `magTheta`::
%	Bending angle of magnet 
% `magLength`::
%	Length of magnet 
%
% OPTIONS
% -------
% The following options can be given as name/value-pairs:
%
% `'refPoint'` (`'center'`)::
%	Reference point where errors are defined. By default the reference point is defined as the center
%   of a straight line connecting the magnet entrance and exit of the ideal beam trajectory. 
%   Optional the reference point can be specified as 'entrance', thus the ideal entrance point of
%   the magnet.
%
% RETURN VALUE
% ------------
% `T1` ::
%	[1x6] array of entrance translation vector to use in lattice element.
% `T2` ::
%	[1x6] array of exit translation vector to use in lattice element.
% `R1` ::
%	[6x6] array of entrance rotation matrix to use in lattice element.
% `R2` ::
%	[6x6] array of exit rotation matrix to use in lattice element.

	% Parse optional arguments
	p = inputParser;
	addOptional(p,'refPoint','center');
	parse(p,varargin{:});
	par = p.Results;
	
	% Define offset vector 
	d0Vector = [dx, dy, dz]';
	
	% Original and x,y,z-axis
	xAxis = [1, 0, 0]';
	yAxis = [0, 1, 0]';
	zAxis = [0, 0, 1]';

	% Radius of curvature
	Rc = magLength/magTheta;

	% Error rotation defined in the error [x/y/z]-frame
	R0 = [ cos(ay)*cos(az)                          ,           -cos(ay)*sin(az)               ,       sin(ay)    ;...
		   cos(az)*sin(ax)*sin(ay) + cos(ax)*sin(az), cos(ax)*cos(az) - sin(ax)*sin(ay)*sin(az), -cos(ay)*sin(ax) ;...
		  -cos(ax)*cos(az)*sin(ay) + sin(ax)*sin(az), cos(az)*sin(ax) + cos(ax)*sin(ay)*sin(az),  cos(ax)*cos(ay) ];

	% Check if errors are given at center of magnet
	switch par.refPoint
		case 'center'
			% Rotation from bending angle
			RB2 = [ cos(magTheta/2), 0 , -sin(magTheta/2) ;...
					0              , 1 ,        0       ;...
					sin(magTheta/2), 0 ,  cos(magTheta/2) ];

			% Rotational-error matrix in the original xyz coordinate system
			RX = RB2 * R0 * RB2';

			if magTheta==0
				OO0 =  (magLength/2) * zAxis ;
				P0P = -(magLength/2) * RX * zAxis;
			else
				OO0 =  Rc*sin(magTheta/2) * RB2 * zAxis;
				P0P = -Rc*sin(magTheta/2) * RX * RB2 * zAxis;
			end

			% Transform offset to magnet entrance
			OP = OO0 + P0P +  RB2 * d0Vector;
		case 'entrance'
			% Keep original error rotation matrix
			RX = R0;
			% Keep original offset vector
			OP = d0Vector;
		otherwise
			error('Unsupported reference point. Allowed are ''center'' or ''entrance''.')
	end
	  
	% Loop over entrance & exit
	for face=1:2
		if face==1
			R = RX;
			% XYZ - axes unit - vectors expressed in the xyz coordinate system
			XaxiSxyz = R * xAxis;
			YaxiSxyz = R * yAxis;
			ZaxiSxyz = R * zAxis;

			LD = ZaxiSxyz' * OP; 
			tmp = OP;

		else
			% Rotation from bending angle
			RB = [ cos(magTheta), 0 , -sin(magTheta) ;...
			                  0 , 1 , 0            ;...
			       sin(magTheta), 0 , cos(magTheta)  ];

			% Rotational - error matrix in the x' y' z' coordinate system
			R = RB' * RX' * RB;
 
			% XYZ - axes unit - vectors expressed in the xyz coordinate system
			XaxiSxyz = RB * xAxis;
			YaxiSxyz = RB * yAxis;
			ZaxiSxyz = RB * zAxis;
			
			if magTheta==0
				OPp = [0, 0, magLength]';
			else
				OPp = [Rc*(cos(magTheta) - 1), 0, magLength*sin(magTheta)/magTheta]';
			end
			OOp = RX * OPp + OP;
			OpPp = (OPp - OOp);
			LD = ZaxiSxyz' * OpPp; 
			tmp = OpPp;
		end


		% 1st translation vector
		tD0 = [-tmp' * XaxiSxyz, 0, -tmp' * YaxiSxyz, 0, 0, 0]';

		% 2nd translation vector 
		T0 = [LD*R(3,1)/R(3,3), R(3,1), LD*R(3,2)/R(3,3), R(3,2),  0 ,LD/R(3,3)]';
		
		% Final translation vector
		T = T0 + tD0;
		
		% Final rotation matrix
		LinMat = [ R(2,2)/R(3,3) , LD*R(2,2)/R(3,3)^2  , -R(1,2)/R(3,3) , -LD*R(1,2)/R(3,3)^2 , 0      , 0 ;...
			          0          , R(1,1)              ,  0             ,  R(2,1)             , R(3,1) , 0 ;...
			      -R(2,1)/R(3,3) , -LD*R(2,1)/R(3,3)^2 ,  R(1,1)/R(3,3) ,  LD*R(1,1)/R(3,3)^2 , 0 , 0      ;...
			          0          , R(1,2)              ,  0             ,  R(2,2)             , R(3,2) , 0 ;...
			          0          , 0                   , 0              , 0                   , 1 , 0      ;
			      -R(1,3)/R(3,3) , -LD*R(1,3)/R(3,3)^2 , -R(2,3)/R(3,3) , -LD*R(2,3)/R(3,3)^2 , 0 , 1      ];
	
		% Asign output
		if face==1
			R1 = LinMat;
			T1 = T;
		else
			R2 = LinMat;
			T2 = T;		
		end
	end
end
