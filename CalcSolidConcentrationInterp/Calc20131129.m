fp=fopen('result.txt','w');
for i=0.0:0.01:1.2
    for j=0.0:0.01:1.2
        %Define the parameter of square
        l=0.1;%length of edge of square
        xmin=i-l/2;
        xmax=i+l/2;
        ymin=j-l/2;
        ymax=j+l/2;
        
        rlb=sqrt(xmin*xmin+ymin*ymin);
        rrb=sqrt(xmax*xmax+ymin*ymin);
        rlt=sqrt(xmin*xmin+ymax*ymax);
        rrt=sqrt(xmax*xmax+ymax*ymax);
        if (rlb<1&&rrb<1&&rlt<1&&rrt<1)
            es=1.0;
        elseif (rlb>1&&rrb>1&&rlt>1&&rrt>1)
            es=0.0;
        else
            %Define grid
            mesh_density=0.001;
            x=xmin:mesh_density:xmax;%-mesh_range:mesh_density:mesh_range;
            y=ymin:mesh_density:ymax;
            [X,Y] = meshgrid(x,y);

            %Define a circle
            r=1;
            oa=0;
            ob=0;
            t=0:0.01:2*pi;
            a=r*sin(t)+oa;
            b=r*cos(t)+ob;
            %plot(a,b)
            in_circle = inpolygon(X,Y,a,b);

            x_polygon = [xmin xmax xmax xmin xmin];
            y_polygon = [ymin ymin ymax ymax ymin];
            %hold on
            %plot(x_polygon,y_polygon);
            %hold off
            in_polygon = inpolygon(X,Y,x_polygon,y_polygon);
            in_both=in_circle & in_polygon;

            %figure 
            %imshow(in_both)

            %figure 
            %imshow(in_polygon)
            %figure 
            %imshow(in_circle)

            Area=sum(sum(in_both))*mesh_density^2;
            es=Area/l/l;
        end
        if (es>1.0)
            es=1.0;
        end
        
        fprintf(1,'%f\t%f\t%f\n',i,j,es);
        fprintf(fp,'%f\t%f\t%f\n',i,j,es);
    end
end
fcolse(fp);