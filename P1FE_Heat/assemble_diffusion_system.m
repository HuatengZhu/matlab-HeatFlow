% Assemble the global diffusion matrix
function DiffMat = assemble_diffusion_system(cell_v,ncell,nvert,vertex);
%% Initialise vectors for sparse matrix
% Evaluate number of non-zeros entries (cf below how many times we do "pos=pos+1")
nz=9*ncell;
IA=zeros(nz,1);
JA=zeros(nz,1);
VA=zeros(nz,1);

% "pos"=position inside the vectors IA, JA, VA that store the entries of A
pos=0;

for i=1:ncell
    % We are collecting vertices for each cell in the following to compute local stiffness matrix
  vertices=[vertex(cell_v{i}(1),1) vertex(cell_v{i}(2),1) vertex(cell_v{i}(3),1);
            vertex(cell_v{i}(1),2) vertex(cell_v{i}(2),2) vertex(cell_v{i}(3),2)]; % this is a 2 X 3 matrix whose j-th column is the coordinate of vertex cell_v{i}(j).
  Sloc = stima(vertices,i);% local stiffness matrix
  
  % Loop over vertices
  for jj=1:3
    jvert = cell_v{i}(jj);
      for kk=1:3
        kvert = cell_v{i}(kk);
        pos=pos+1;
        IA(pos) = jvert;
        JA(pos) = kvert;
        VA(pos) = Sloc(jj,kk);%It does not overwrite by definition
      end
    %end
  end
end

%% Creation of the sparse matrix
DiffMat=sparse(IA(1:pos),JA(1:pos),VA(1:pos),nvert,nvert);

