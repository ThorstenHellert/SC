function [T1,T2,R1,R2] = SCgetTransformation(dx,dy,dz,ax,ay,az,magTheta,magLength)
% SCgetTransformation
% ===================
%
% NAME
% ----
% SCgetTransformation - calculates the translation and rotation elements for given misalignments
%
% SYNOPSIS
% --------
% `[T1,T2,R1,R2] = SCgetTransformation(dx,dy,dz,ax,ay,az,thetaB,lenB)`
%
%
% DESCRIPTION
% -----------
% This function calculates AT's linear rotation and translation elements `T1`, `T2`, and
% `R1`, `R2` for magnets with aribtrary offsets and angles, see Tech-Note XXX.
%
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

	
	
	% Define offset vector 
	dVector = [dx, dy, dz]';
	% Original and x,y,z-axis
	xAxis = [1, 0, 0]';
	yAxis = [0, 1, 0]';
	zAxis = [0, 0, 1]';

	% Calculate relevant sin and cos
	cx = cos(ax); 
	sx = sin(ax);
	cy = cos(ay); 
	sy = sin(ay);
	cz = cos(az); 
	sz = sin(az);


	RX = [ cy*cz,           -cy*sz,             sy    ;...
		   cz*sx*sy + cx*sz, cx*cz - sx*sy*sz, -cy*sx ;...
		  -cx*cz*sy + sx*sz, cz*sx + cx*sy*sz,  cx*cy ];

	% Loop over entrance & exit
	for face=1:2
		if face==1
			R = RX;
			Xaxis = R * xAxis;
			Yaxis = R * yAxis;
			Zaxis = R * zAxis;
			Xaxis2 = Xaxis;
			Yaxis2 = Yaxis;
			OP = dVector;
			LD = Zaxis' * dVector; 
			tmp = OP;

		else
			RB = [ cos(magTheta), 0 , -sin(magTheta) ;...
			                0 , 1 , 0            ;...
			       sin(magTheta), 0 , cos(magTheta)  ];

			RBT = RB';
			R = RBT * RX' * RB;
 
			Xaxis2 = RB * xAxis;
			Yaxis2 = RB * yAxis;
			Zaxis2 = RB * zAxis;
			
			if magTheta==0
				OPp = [0, 0, magLength]';
			else
				OPp = [magLength/magTheta*(cos(magTheta) - 1), 0, magLength*sin(magTheta)/magTheta]';
			end
			OOp = RX * OPp + dVector;
			OpPp = (OPp - OOp);
			LD = Zaxis2' * OpPp; 
			tmp = OpPp;
		end


		% Intial translation vector
		tD0 = [-tmp' * Xaxis2, 0, -tmp' * Yaxis2, 0, 0, 0]';

		% Final rotation matrix
		R = [ R(2,2)/R(3,3) , LD*R(2,2)/R(3,3)^2  , -R(1,2)/R(3,3) , -LD*R(1,2)/R(3,3)^2 , 0      , 0 ;...
			     0          , R(1,1)              ,  0             ,  R(2,1)             , R(3,1) , 0 ;...
			 -R(2,1)/R(3,3) , -LD*R(2,1)/R(3,3)^2 ,  R(1,1)/R(3,3) ,  LD*R(1,1)/R(3,3)^2 , 0 , 0      ;...
				 0          , R(1,2)              ,  0             ,  R(2,2)             , R(3,2) , 0 ;...
				 0          , 0                   , 0              , 0                   , 1 , 0      ;
			 -R(1,3)/R(3,3) , -LD*R(1,3)/R(3,3)^2 , -R(2,3)/R(3,3) , -LD*R(2,3)/R(3,3)^2 , 0 , 1      ];

		% Exit translation vector 
		T0 = [LD*R(3, 1)/R(3, 3), R(3, 1), LD*R(3, 2)/R(3, 3), R(3, 2),  0 ,LD/R(3, 3)]';

		% Final translation vector
		T = T0 + tD0;
		
		% Asign output
		if face==1
			R1 = R;
			T1 = T;
		else
			R2 = R;
			T2 = T;		
		end
	end
end
