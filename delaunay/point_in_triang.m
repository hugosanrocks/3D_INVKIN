function res = poin_in_triang(px,py,p0x,p0y,p1x,p1y,p2x,p2y)

Area = 0.5 *(-p1y*p2x + p0y*(-p1x + p2x) + p0x*(p1y - p2y) + p1x*p2y);

s = 1/(2*Area)*(p0y*p2x - p0x*p2y + (p2y - p0y)*px + (p0x - p2x)*py);
t = 1/(2*Area)*(p0x*p1y - p0y*p1x + (p0y - p1y)*px + (p1x - p0x)*py);

        det2 = px*(p1y-p2y)-py*(p1x-p2x)+(p1x*p2y-p2x*p1y);
        det1 = px*(p0y-p1y)-py*(p0x-p1x)+(p0x*p1y-p1x*p0y);
        det3 = px*(p0y-p2y)-py*(p0x-p2x)+(p0x*p2y-p2x*p0y);

        dp1 = sqrt((p0x-p1x)**2.+(p0y-p1y)**2.);
        dp2 = sqrt((p1x-p2x)**2.+(p1y-p2y)**2.);
        dp3 = sqrt((p0x-p2x)**2.+(p0y-p2y)**2.);

        d1 = sqrt((px-p0x)^2+(py-p0y)^2);
        d2 = sqrt((px-p1x)^2+(py-p1y)^2);
        d3 = sqrt((px-p2x)^2+(py-p2y)^2);

         if ( ( s > 0 ) && ( t > 0. ) && ( (1-s-t ) > 0 ) )
          res = 1;
         else
          res = 0;
          if ( (d1 == 0) || (d2 == 0) || (d3 == 0) )
           res = 1;
          elseif ( ((det1 == 0) && ((d1 < dp1) && (d2 < dp1)) ) || ...
                   ((det2 == 0.) && ((d2 < dp2) && (d3 < dp2)) ) || ...
                   ((det3 == 0.) && ((d3 < dp3) && (d1 < dp3)) ) )
           res = 1;
          end
        end

return
