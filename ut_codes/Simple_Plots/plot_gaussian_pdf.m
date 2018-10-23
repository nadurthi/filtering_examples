function plot_gaussian_pdf(mux,Px,type)
sigs=eig(Px);
[xx,yy]=meshgrid(linspace(mux(1)-3*sigs(1),mux(1)+3*sigs(1),25),linspace(mux(2)-3*sigs(2),mux(2)+3*sigs(2),25));
f=zeros(size(xx));
for i=1:1:size(xx,1)
    for j=1:1:size(xx,2)
        f(i,j)=mvnpdf([xx(i,j),yy(i,j)],mux(:)',Px);
    end
end
if strcmp(type,'contours')
contour(xx,yy,f,20)
end

if strcmp(type,'surface')
mesh(xx,yy,f)
end

end