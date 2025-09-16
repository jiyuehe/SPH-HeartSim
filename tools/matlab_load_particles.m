% to run this code in Matlab:
% in the "EDITOR" tab at top, click the green triangle "Run" button

temp_dir = pwd;
[parent_dir, ~, ~] = fileparts(temp_dir);
data_dir = fullfile(parent_dir, '/build/sim/bin/reload');

file_name = 'HeartModel_rld.xml';
cd(data_dir);
fid = fopen(file_name,'r');

s = fgetl(fid);
while ~contains(s, '<particle OriginalID=')
    s = fgetl(fid);
end

OriginalID = [];
HeartModelPart2ID = [];
VolumetricMeasure = [];
Phi = [];
Position = [];
Fiber = [];
Sheet = [];
while ~contains(s, '</particles>')
    id1 = strfind(s,'OriginalID="') + length('OriginalID="');
    id2 = strfind(s,'" HeartModelPart2ID') - 1;
    OriginalID(end+1) = str2double(s(id1:id2));

    id1 = strfind(s,'HeartModelPart2ID="') + length('HeartModelPart2ID="');
    id2 = strfind(s,'" VolumetricMeasure') - 1;
    HeartModelPart2ID(end+1) = str2double(s(id1:id2));

    id1 = strfind(s,'VolumetricMeasure="') + length('VolumetricMeasure="');
    id2 = strfind(s,'" Phi') - 1;
    VolumetricMeasure(end+1) = str2double(s(id1:id2));

    id1 = strfind(s,'Phi="') + length('Phi="');
    id2 = strfind(s,'" Position') - 1;
    Phi(end+1) = str2double(s(id1:id2));

    id1 = strfind(s,'Position="') + length('Position="');
    id2 = strfind(s,'" Fiber') - 1;
    Position(end+1,:) = str2num(s(id1:id2));

    id1 = strfind(s,'Fiber="') + length('Fiber="');
    id2 = strfind(s,'" Sheet') - 1;
    Fiber(end+1,:) = str2num(s(id1:id2));

    id1 = strfind(s,'Sheet="') + length('Sheet="');
    id2 = strfind(s,'"/>') - 1;
    Sheet(end+1,:) = str2num(s(id1:id2));

    s = fgetl(fid);
end

fclose(fid);

% plot
% figure;
% plot(OriginalID,'.');
% 
% figure;
% plot(HeartModelPart2ID,'.');
% 
% figure;
% plot(VolumetricMeasure,'.');
% 
% figure;
% plot(Phi,'.');

figure;
id1 = HeartModelPart2ID == 0;
id2 = HeartModelPart2ID == 2;
scatter3(Position(id1,1),Position(id1,2),Position(id1,3),'b','.');
hold on;
scatter3(Position(id2,1),Position(id2,2),Position(id2,3),'r','.');
hold off;
axis vis3d equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
rotate3d on;

% figure;
% quiver3(Position(:,1),Position(:,2),Position(:,3),Fiber(:,1),Fiber(:,2),Fiber(:,3));
% axis vis3d equal;
% rotate3d on;
% 
% figure;
% quiver3(Position(:,1),Position(:,2),Position(:,3),Sheet(:,1),Sheet(:,2),Sheet(:,3));
% axis vis3d equal;
% rotate3d on;
% 
% dot_product = zeros(size(Fiber,1),1);
% for n = 1:size(Fiber,1)
%     dot_product(n) = dot(Fiber(1,:),Sheet(1,:));
% end
% figure;
% plot(dot_product,'.'); % Sheet is perpendicular to Fiber
