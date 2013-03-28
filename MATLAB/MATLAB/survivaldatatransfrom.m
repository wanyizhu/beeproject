   for i=1:69
for j=1:19
x(i+1,j)=death(j,)-death(j,i);
end


for i=1:19
x(i+1,:)=condeath(i+1,:)-condeath(i,:);
end

for i=1:69
for j=1:20
t{i,j}=repmat(days(i,j),dt(i,j),1);
mat(i,:)=vertcat(t{i,:});
end 
end 



st=death(1,:)
for i=1:20
end
B = repmat(st,20,69)
B = repmat(st,20,1)
death
condeath=B-death
for i=1:19
ttt(i+1,:)=condeath(i+1,:)-condeath(i,:);
end
ttt
dt=ttt';
for i=1:69
for j=1:20
t{i,j}=repmat(days(j),dt(i,j),1);
end
end
t
vertcat(t{:,1})

for i=1:29
for j=1:6
t{i,j}=repmat(days(j),dt(i,j),1);
end
end

a1=vertcat(t{1,:});a2=vertcat(t{2,:});a3=vertcat(t{3,:});a4=vertcat(t{4,:});a5=vertcat(t{5,:});a6=vertcat(t{6,:});a7=vertcat(t{7,:});a8=vertcat(t{8,:});a9=vertcat(t{9,:});a10=vertcat(t{10,:});