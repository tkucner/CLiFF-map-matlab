function [] = SaveXML(obj,FileName)

docNode = com.mathworks.xml.XMLUtils.createDocument('map');
map = docNode.getDocumentElement;
map.setAttribute('version','0.1');

parameters = docNode.createElement('parameters');

x_min = docNode.createElement('x_min');
x_min.appendChild(docNode.createTextNode(num2str(obj.GridParameters(1))));
parameters.appendChild(x_min);

x_max = docNode.createElement('x_max');
x_max.appendChild(docNode.createTextNode(num2str(obj.GridParameters(2))));
parameters.appendChild(x_max);

y_min = docNode.createElement('y_min');
y_min.appendChild(docNode.createTextNode(num2str(obj.GridParameters(3))));
parameters.appendChild(y_min);

y_max = docNode.createElement('y_max');
y_max.appendChild(docNode.createTextNode(num2str(obj.GridParameters(4))));
parameters.appendChild(y_max);

step = docNode.createElement('step');
step.appendChild(docNode.createTextNode(num2str(obj.GridParameters(5))));
parameters.appendChild(step);

radious = docNode.createElement('radious');
radious.appendChild(docNode.createTextNode(num2str(obj.Radious)));
parameters.appendChild(radious);

wind = docNode.createElement('wind');
wind.appendChild(docNode.createTextNode(num2str(obj.Wind)));
parameters.appendChild(wind);

map.appendChild(parameters);

locations = docNode.createElement('locations');

for i=1:numel(obj.Batches)
    
    location = docNode.createElement('location');
    
    id = docNode.createElement('id');
    id.appendChild(docNode.createTextNode(num2str(obj.Batches(i).ID)));
    location.appendChild(id);
    
    pose = docNode.createElement('pose');
    
    x_pose = docNode.createElement('x');
    x_pose.appendChild(docNode.createTextNode(num2str(obj.Batches(i).Pose(1))));
    pose.appendChild(x_pose);
    
    y_pose = docNode.createElement('y');
    y_pose.appendChild(docNode.createTextNode(num2str(obj.Batches(i).Pose(2))));
    pose.appendChild(y_pose);
    
    location.appendChild(pose);
    
    for j=1:numel(obj.Batches(i).P)
        
        distribution = docNode.createElement('distribution');
        
        Cluster = docNode.createElement('Cluster');
        Cluster.appendChild(docNode.createTextNode(num2str(obj.Batches(i).MapCluster_Direction(j))));
        distribution.appendChild(Cluster);        
        
        P = docNode.createElement('P');
        P.appendChild(docNode.createTextNode(num2str(obj.Batches(i).P(j))));
        distribution.appendChild(P);
        
        M = docNode.createElement('M');
        
        th = docNode.createElement('th');
        th.appendChild(docNode.createTextNode(num2str(obj.Batches(i).Mean(j,1))));
        M.appendChild(th);
        
        r = docNode.createElement('r');
        r.appendChild(docNode.createTextNode(num2str(obj.Batches(i).Mean(j,2))));
        M.appendChild(r);
        
        distribution.appendChild(M);
        
        Cov = docNode.createElement('Cov');
        e_11 = docNode.createElement('e_11');
        e_11.appendChild(docNode.createTextNode(num2str(obj.Batches(i).Cov(1,1,j))));
        Cov.appendChild(e_11);
        
        e_12 = docNode.createElement('e_12');
        e_12.appendChild(docNode.createTextNode(num2str(obj.Batches(i).Cov(1,2,j))));
        Cov.appendChild(e_12);
        
        e_21 = docNode.createElement('e_21');
        e_21.appendChild(docNode.createTextNode(num2str(obj.Batches(i).Cov(2,1,j))));
        Cov.appendChild(e_21);
        
        e_22 = docNode.createElement('e_22');
        e_22.appendChild(docNode.createTextNode(num2str(obj.Batches(i).Cov(2,2,j))));
        Cov.appendChild(e_22);
        
        distribution.appendChild(Cov);
        
        location.appendChild(distribution);        
    end
    locations.appendChild(location);
end
map.appendChild(locations);

locations_sparse = docNode.createElement('locations_sparse');

for i=1:numel(obj.BatchesSparse)
    
    location = docNode.createElement('location');
    
    id = docNode.createElement('id');
    id.appendChild(docNode.createTextNode(num2str(obj.BatchesSparse(i).ID)));
    location.appendChild(id);
    
    pose = docNode.createElement('pose');
    
    x_pose = docNode.createElement('x');
    x_pose.appendChild(docNode.createTextNode(num2str(obj.BatchesSparse(i).Pose(1))));
    pose.appendChild(x_pose);
    
    y_pose = docNode.createElement('y');
    y_pose.appendChild(docNode.createTextNode(num2str(obj.BatchesSparse(i).Pose(2))));
    pose.appendChild(y_pose);
    
    location.appendChild(pose);
    
    for j=1:numel(obj.BatchesSparse(i).P)
        
        distribution = docNode.createElement('distribution');
        
        Cluster = docNode.createElement('Cluster');
        Cluster.appendChild(docNode.createTextNode(num2str(obj.BatchesSparse(i).MapCluster_Direction(j))));
        distribution.appendChild(Cluster);        
        
        P = docNode.createElement('P');
        P.appendChild(docNode.createTextNode(num2str(obj.BatchesSparse(i).P(j))));
        distribution.appendChild(P);
        
        M = docNode.createElement('M');
        
        th = docNode.createElement('th');
        th.appendChild(docNode.createTextNode(num2str(obj.BatchesSparse(i).Mean(j,1))));
        M.appendChild(th);
        
        r = docNode.createElement('r');
        r.appendChild(docNode.createTextNode(num2str(obj.BatchesSparse(i).Mean(j,2))));
        M.appendChild(r);
        
        distribution.appendChild(M);
        
        Cov = docNode.createElement('Cov');
        e_11 = docNode.createElement('e_11');
        e_11.appendChild(docNode.createTextNode(num2str(obj.BatchesSparse(i).Cov(1,1,j))));
        Cov.appendChild(e_11);
        
        e_12 = docNode.createElement('e_12');
        e_12.appendChild(docNode.createTextNode(num2str(obj.BatchesSparse(i).Cov(1,2,j))));
        Cov.appendChild(e_12);
        
        e_21 = docNode.createElement('e_21');
        e_21.appendChild(docNode.createTextNode(num2str(obj.BatchesSparse(i).Cov(2,1,j))));
        Cov.appendChild(e_21);
        
        e_22 = docNode.createElement('e_22');
        e_22.appendChild(docNode.createTextNode(num2str(obj.BatchesSparse(i).Cov(2,2,j))));
        Cov.appendChild(e_22);
        
        distribution.appendChild(Cov);
        
        location.appendChild(distribution);        
    end

    locations_sparse.appendChild(location);    
end
map.appendChild(locations_sparse);

xmlwrite(FileName,docNode);
end