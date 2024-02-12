clc
clear all
close all
%calculate velocity script
concentrations = [{'Dead'},{'RES'},{'Water'},{'0.0'},{'0.1'},{'1.0'},{'10'},{'100'}];

for i = 1:length(concentrations)
    index = 1;
    conc_dir = ['./',concentrations{i}];
    files = dir(conc_dir);
    disp(concentrations{i});
    excel_files = cell(1);
    for j = 3:length(files)
        inst_dir = [conc_dir,'/',files(j).name];
        sub_files = dir(inst_dir);
        for z = 3:length(sub_files)
            [path,name,ext] = fileparts([inst_dir,'/',sub_files(z).name]);
            if(strcmp(ext,'.xls'))%&&strcmp('09_16',name(1:5)));
                excel_files(index) = {[path,'/',name,ext]};
                disp([name,ext]);
                index = index + 1;
            end
        end
    end
    
    if(~isempty(excel_files{1}))
        total_V = zeros(599,13);
        date_V = total_V;
        num_files = 0;
        date_files = 0;
        for j = 1:length(excel_files)
            [~,name,~] = fileparts(excel_files{j});
            data = xlsread(excel_files{j});
            X = data(1:2:end,:);
            Y = data(2:2:end,:);
            V = ((diff(X,1,1)).^2+(diff(Y,1,1)).^2).^(.5);
            %multiplying to convert to um
            V = V.*(400/150);
            K = menger_curve(data);
            xlswrite(excel_files{j},[transpose(.1:.1:length(V)/10),V],'Velocity');
            xlswrite(excel_files{j},[transpose(.1:.1:length(K)/10),K],'Curvature');
            if(length(V)==length(total_V))
                total_V = total_V+V;
                date_V = date_V+V;
                date_files = date_files+1;
                num_files = num_files+1;
            else
                disp(excel_files{j})
            end
            
            if(j<length(excel_files))
                [~,next_name,~] = fileparts(excel_files{j+1});
                if(~strcmp(name(1:5),next_name(1:5)))
                    date_V = date_V./date_files;
                    t = transpose(.1:.1:599/10);
                    xlswrite([conc_dir,'/',concentrations{i},'_',name(1:5),'_averaged.xls'],[t,date_V],'Velocity');
                    xlswrite([conc_dir,'/',concentrations{i},'_',name(1:5),'_averaged.xls'],[{'Number of worms'},{date_files}],'n');
                    date_V = zeros(599,13);
                    date_files = 0;
                end
            else
                date_V = date_V./date_files;
                t = transpose(.1:.1:599/10);
                xlswrite([conc_dir,'/',concentrations{i},'_',name(1:5),'_averaged.xls'],[t,date_V],'Velocity');
                xlswrite([conc_dir,'/',concentrations{i},'_',name(1:5),'_averaged.xls'],[{'Number of worms'},{date_files}],'n');
                date_V = zeros(599,13);
                date_files = 0;
            end
        end
        total_V = total_V./num_files;
        t = transpose(.1:.1:599/10);

        xlswrite([conc_dir,'/',concentrations{i},'_averaged.xls'],[t,total_V],'Velocity');
        xlswrite([conc_dir,'/',concentrations{i},'_averaged.xls'],[{'Number of worms'},{num_files}],'n');
    end
end

%%
sections = [{'Head'},{'Mid'},{'Tail'}];
indxs = [1,7,13];
combined_head = zeros(599,length(concentrations));
combined_head = [transpose(.1:.1:599/10),combined_head];
combined_mid = combined_head;
combined_tail = combined_head;
for i = 1:length(concentrations)
    data = xlsread(['./',concentrations{i},'/',concentrations{i},'_averaged.xls'],'Velocity');
    if(size(data,2)>1)
        combined_head(:,i+1) = data(:,2);
        combined_mid(:,i+1) = data(:,8);
        combined_tail(:,i+1) = data(:,14);
    end
end
combined_head = [[sections(1),concentrations];num2cell(combined_head)];
combined_mid = [[sections(2),concentrations];num2cell(combined_mid)];
combined_tail = [[sections(3),concentrations];num2cell(combined_tail)];

xlswrite('Combined.xls',combined_head,sections{1});
xlswrite('Combined.xls',combined_mid,sections{2});
xlswrite('Combined.xls',combined_tail,sections{3});
