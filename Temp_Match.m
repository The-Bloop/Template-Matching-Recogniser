classdef Temp_Match
    properties
       img;
       ind;
       temp_img;
       temp_ind;
       tp;
       tn;
       fp;
       fn;
       hd_score;
       mhd_score;
       tani_score;
       yule_score;
       score;
    end
    
    methods
        
        function SetTemplate(this,tempstr)
            this.temp_img = imread(tempstr);
%             f2 = figure; figure(f2); imshow(obj.temp_img,[]);
            a = find(this.temp_img == 0);
            [row,col] = ind2sub([48,48],a);
            this.temp_ind = [row,col];
        end
        
        function score = GetScore(this)
           score = this.score; 
        end
        
        function img = Convert2Img(this, points)
            X = points(:,1);	X = X-min(X);
            Y = points(:,2);	Y = Y-min(Y);
            plot(X,Y,'o');
            
%             r = sqrt(X.^2 + y.^2)
%             theta = atan2(Y,X)
%             plot(theta,r)
            
            %rx = diff(X);
            %ry = diff(Y);
            %r =cumsum(sqrt(rx.^2 + ry.^2));
            %theta = atan([rx,ry]);
            %f1 = figure;
            %figure(f1);
            %plot(theta,r,'*');
            
            H = max(Y);
            L = max(X);
            d = abs(H-L)/2;
            if H>L
                X2 = (X+d/2)/(L+d)*47;
                Y2 = Y/H*47;
            else
                Y2 = (Y+d/2)/(H+d)*47;
                X2 = X/L*47;
            end
            X3 = round(X2);
            Y3 = round(Y2);
            
            
            %t_max = max(abs(theta))
            %theta = theta + t_max
            
            X4 = abs(X3) + 1;
            Y4 = abs(Y3 - 47) + 1;
            
            p = [X4 Y4];
            q = unique(p,'rows');
            X5 = q(:,1);
            Y5 = q(:,2);
            
            
            
            idx1 = sub2ind([48 48], Y4, X4);
            idx1 = unique(idx1);
            idx1 = sort(idx1);
            idx2 = sub2ind([48 48], Y5, X5);
            idx2 = sort(idx2);
            m = mean(idx1==idx2);
            img = ones(48) .* 255;
            img(idx2) = 0;
            this.img = img;
            a = find(img == 0);
            [row,col] = ind2sub([48,48],a);
            this.ind = [row,col];
%             f2 = figure;
%             figure(f2);
%             imshow(img,[]);
        end
        
        function score = HD(this, temp, img)
            temp_img = temp;
%             f2 = figure; figure(f2); imshow(obj.temp_img,[]);
            b = find(temp_img == 0);
            [row,col] = ind2sub([48,48],b);
            temp_ind = [row,col];
            
            a = find(img == 0);
            [row,col] = ind2sub([48,48],a);
            ind = [row,col];
            
            A = ind;
            B = temp_ind;
            len = length(A);
            min_array = zeros(len,1);
            for i=1:length(A)
                d = A(i,:) - B;
                d = sqrt(d(:,1).^2 + d(:,2).^2);
                dmin = min(d);
                min_array(i) = dmin;
            end
            u_min_array = unique(min_array);
            k = round(length(u_min_array)*(6/100));
            new_len = length(u_min_array) - k;
            u_min_array = sort(u_min_array);
            new_min_array = u_min_array(1:new_len);
            score = max(new_min_array);   
        end
        
        function score = MHD(this, temp, img)
            A = img;
            B = temp;
            scoreAB = this.mhd(A,B);
            scoreBA = this.mhd(B,A);
            score_array = [scoreAB, scoreBA];
            score = max(score_array);
        end
        
        function [scoret,scorey] = TnY(this, temp, img)
            temp_img = temp;
            b = find(temp_img == 0);
            [row,col] = ind2sub([48,48],b);
            temp_ind = [row,col];
            
            a = find(img == 0);
            [row,col] = ind2sub([48,48],a);
            ind = [row,col];
            
            ca = find(img > 0);
            
            [TP,TN,FP,FN] = get_data(this, temp_img, img);
            
            T_ab = TP/(TP + FN + FP);
            Tc_ab = FN/(FN + FP + FN);
            
            alpha = 0.75 - 0.25*(0.25*((2*TP) + FP + FN))/(2*48*48);
            
            scoret = (alpha * T_ab) + ((1-alpha)*Tc_ab);
            
            y_n = (TP*TN) - (FN*FP);
            
            y_d = (TP*TN) + (FN*FP);
            
            scorey = y_n/y_d;           
        end
        
        function [TP,TN,FP,FN] = get_data(this, temp, img)
            temp_img = temp;
            b = find(temp_img == 0);
            [row,col] = ind2sub([48,48],b);
            temp_ind = [row,col];
            
            a = find(img == 0);
            [row,col] = ind2sub([48,48],a);
            ind = [row,col];
            
            x = size(ind);
            y = size(b);
            TP = this.get_intst(ind,b);
            FN = length(a) - TP;
            FP = length(b) - TP;
            TN = (48*48) - TP - FN - FP;   
        end
        
        function count = get_intst(this,a,b)
           len = length(a);
           count = 0;
           for i=1:len
               na1 = ((a(i,1)-4):1:(a(i,1)+4))';
               na2 = ((a(i,2)-4):1:(a(i,2)+4))';
               na = na1 * na2';
               na = reshape(na, [81,1]);
               [val,index] = intersect(b, na);
               count = count + length(val);
               for j=1:length(val)
                   m = find(b == val(j));
                   b(m) = [];
               end
           end
        end
        
        function score = mhd(this, temp, img)
            temp_img = temp;
%             f2 = figure; figure(f2); imshow(obj.temp_img,[]);
            a = find(temp_img == 0);
            [row,col] = ind2sub([48,48],a);
            temp_ind = [row,col];
            
            a = find(img == 0);
            [row,col] = ind2sub([48,48],a);
            ind = [row,col];
            
            A = ind;
            B = temp_ind;
            len = length(A);
            min_array = zeros(len,1);
            for i=1:length(A)
                d = A(i,:) - B;
                d = sqrt(d(:,1).^2 + d(:,2).^2);
                dmin = min(d);
                min_array(i) = dmin;
            end
            score = mean(min_array);
        end
        
        
    end 
end
