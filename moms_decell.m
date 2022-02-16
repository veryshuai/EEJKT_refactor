function [lambda_f, lambda_h, c_val_f, c_val_h,value_f,value_h]  = moms_decell(lambda_f, lambda_h, c_val_f, c_val_h,value_f,value_h)
% This function takes policy cells and makes them into multi-dimensional arrays

    %% Decell array 
    lambda_f_temp = zeros(size(lambda_f,1),size(lambda_f,2),size(lambda_f,3),size(lambda_f,4),size(lambda_f{1,1,1,1},1),size(lambda_f{1,1,1,1},2));
    lambda_h_temp = zeros(size(lambda_h,1),size(lambda_h,2),size(lambda_h,3),size(lambda_h{1,1,1},1),size(lambda_h{1,1,1},2));
    value_f_temp = zeros(size(value_f,1),size(value_f,2),size(value_f,3),size(value_f,4),size(value_f{1,1,1,1},1),size(value_f{1,1,1,1},2));
    value_h_temp = zeros(size(value_h,1),size(value_h,2),size(value_h,3),size(value_h{1,1,1},1),size(value_h{1,1,1},2));
    c_val_h_temp = zeros(size(c_val_h,1),size(c_val_h{1},1),size(c_val_h{1},2));
    c_val_f_temp = zeros(size(c_val_f,1),size(c_val_f{1},1),size(c_val_f{1},2));
    for m = 1:size(lambda_f,1)
        for n = 1:size(lambda_f,2)
            for o = 1:size(lambda_f,3)
                for p = 1:size(lambda_f,4)
                    lambda_f_temp(m,n,o,p,:,:) = lambda_f{m,n,o,p}(:,:);
                    value_f_temp(m,n,o,p,:,:) = value_f{m,n,o,p}(:,:);
                end
            end
        end
    end
    for m = 1:size(lambda_h,1)
        for n = 1:size(lambda_h,2)
            for o = 1:size(lambda_h,3)
                lambda_h_temp(m,n,o,:,:) = lambda_h{m,n,o}(:,:);
                value_h_temp(m,n,o,:,:) = value_h{m,n,o}(:,:);
            end
        end
    end
    for m = 1:size(c_val_h,1)
                c_val_h_temp(m,:,:) = c_val_h{m}(:,:);
    end
    for m = 1:size(c_val_f,1)
                c_val_f_temp(m,:,:) = c_val_f{m}(:,:);
    end
    lambda_f = lambda_f_temp; %(succ,trial,common succ rate (defunct), network size, prod of firm, macro shock) 
    lambda_h = lambda_h_temp; %(common succ rate (defunct), known succ rate, network size)
    value_f = value_f_temp; %(succ,trial,common succ rate (defunct), network size, prod of firm, macro shock) 
    value_h = value_h_temp; %(common succ rate (defunct), known succ rate, network size)
    c_val_h = c_val_h_temp; %(demand shock, prod of firm, macro shock)
    c_val_f = c_val_f_temp; %(demand shock, prod of firm, macro shock)
