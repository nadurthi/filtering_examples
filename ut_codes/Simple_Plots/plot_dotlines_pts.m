function plot_dotlines_pts(x,w,col)
if size(x,2)==1
  for i=1:1:size(x,1)
      plot([x(i),x(i)],[0,w(i)],col,'Linewidth',1.5) 
  end
  return
end

for i=1:1:size(x,1)
    plot3([x(i,1),x(i,1)],[x(i,2),x(i,2)],[0,w(i)],col,'Linewidth',1.5);
end

