% Returns Beckmann & haidvogel; Haney numbers in that order
%           [r_x0,r_x1] = stiffness(h,zrmat)

function [r_x0,r_x1] = stiffness(h,zrmat)

    % Grid metrics - From Utility/stiffness.F
    % beckmann & haidvogel (1993) number 
    r_x0_X = abs(diff(h,1,1)./avg1(h,1))/2;
    r_x0_Y = abs(diff(h,1,2)./avg1(h,2))/2;
    r_x0 = max([max(r_x0_X(:)) max(r_x0_Y(:))]);

    clear r_x0_X r_x0_y

    % haney (1991) number
    r_x1_x = nan(size(zrmat));
    r_x1_y = nan(size(zrmat));

    for k=2:size(zrmat,3)
          for i=2:size(zrmat,1)
            r_x1_x(i,:,k) = (zrmat(i,:,k) - zrmat(i-1,:,k) + zrmat(i,:,k-1) - zrmat(i-1,:,k-1)) ...
                            ./ (zrmat(i,:,k) + zrmat(i-1,:,k) - zrmat(i,:,k-1) - zrmat(i-1,:,k-1));
          end

          for j=2:size(zrmat,2)
            r_x1_y(:,j,k) = (zrmat(:,j,k) - zrmat(:,j-1,k) + zrmat(:,j,k-1) - zrmat(:,j-1,k-1)) ...
                            ./ (zrmat(:,j,k) + zrmat(:,j-1,k) - zrmat(:,j,k-1) - zrmat(:,j-1,k-1));
          end
    end
    r_x1 = nanmax(abs( [r_x1_x(:);r_x1_y(:)] ));
    clear r_x1_x r_x1_y