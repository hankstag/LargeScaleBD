%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code implementing the paper "Large-Scale Bounded Distortion Mappings".
% Disclaimer: The code is provided as-is for academic use only and without any guarantees. 
%             Please contact the author to report any bugs.
% Written by Shahar Kovalsky (http://www.wisdom.weizmann.ac.il/~shaharko/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bd_map = bd_map(filename)
%% init
rng(1)
%addpath('mex');


%% parameters
n = 100; % problem size
sigma = .2; % noise level (for initial map)
K = 1e10; % conformal distortion bound
lb = -1; % lower bound on SVs (-1 = disabled)
ub = -1; % upper bound on SVs (-1 = disabled)
iter_max = 1000; % maximal number of BD projection iterations
tol_err = 1e-15; % tolerance for stopping BD projection iterations
use_weighted_metric = false; % use a weighted metric?


%% generate problem
% generate a noisy surface
filename
obj = readObj("models/"+filename);
 
V = obj.v;
F = obj.f.v;

% some constants
dim = size(F,2)-1;
n_vert = size(V,1);
n_tri = size(F,1);

% initial map
x0 = obj.vt;

% setup linear constraints (fix centroid)
n = 6;
s = size(x0,1);
idx = strfind(filename,".obj");
model_name = extractBefore(filename,idx);
eq_l = csvread("constraints/"+model_name+"_pi.csv");

eq_lhs = zeros(n,s*dim);
a = eq_l(1,1);
eq_lhs(1,a) = 1;
eq_lhs(2,a+s) = 1;
b = eq_l(2,1);
eq_lhs(3,b) = 1;
eq_lhs(4,b+s) = 1;
c = eq_l(3,1);
eq_lhs(5,c) = 1;
eq_lhs(6,c+s) = 1;
eq_rhs = zeros(n,1);
eq_r = csvread("constraints/"+model_name+"_p.csv");
eq_rhs = [eq_r(1,:)';eq_r(2,:)';eq_r(3,:)']
% eq_rhs(1) = 1;
% eq_rhs(2) = 1;
% eq_rhs(3) = -1;
% eq_rhs(4) = 1;
% eq_rhs(5) = 0;
% eq_rhs(6) = -1;
%eq_rhs(10,1) = -1;
% eq_lhs = kron(eye(dim),ones(1,n_vert))/n_vert;
% eq_rhs = eq_lhs*colStack(x0);

% plot surface
% figure;
% patch('vertices',V,'faces',F,'facecolor','w');
% axis equal; axis off;
% cameratoolbar;
% cameratoolbar('SetCoordSys','none');
% title('surface');


%% solve problem
% setup BD solver
solver_bd = SolverProjectorBD(F, V, eq_lhs, eq_rhs, K, lb, ub, x0, SolverProjectorModeEnum.Tangent, use_weighted_metric);

% plot initial map
% figure;
% solver_bd.visualize();
% title('Initial Map');
% hA(1) = gca;

% run solver
solver_bd.solve(iter_max, tol_err); % solve BD projection

%%
OBJ.vertices = V;
OBJ.vertices_texture = solver_bd.y;
OBJ.objects(1).type='f';
OBJ.objects(1).data.vertices=F;
OBJ.objects(1).data.texture =F;
OBJ.objects(1).data.normal = F;
write_wobj(OBJ,"LGBD/"+filename);
%
% OBJ struct containing:
%
% OBJ.vertices : Vertices coordinates
% OBJ.vertices_texture: Texture coordinates 
% OBJ.vertices_normal : Normal vectors
% OBJ.vertices_point  : Vertice data used for points and lines   
% OBJ.material : Parameters from external .MTL file, will contain parameters like
%           newmtl, Ka, Kd, Ks, illum, Ns, map_Ka, map_Kd, map_Ks,
%           example of an entry from the material object:
%       OBJ.material(i).type = newmtl
%       OBJ.material(i).data = 'vase_tex'
% OBJ.objects  : Cell object with all objects in the OBJ file, 
%           example of a mesh object:
%       OBJ.objects(i).type='f'               
%       OBJ.objects(i).data.vertices: [n x 3 double]
%       OBJ.objects(i).data.texture:  [n x 3 double]
%       OBJ.objects(i).data.normal:   [n x 3 double]


% disp(solver_bd.y(a,:));
% disp(solver_bd.y(b,:));
% disp(solver_bd.y(c,:));
% % plot output map
% figure;
% solver_bd.visualize();
% title('Output Map');
% hA(2) = gca;
% linkaxes(hA);
%disp(a)
