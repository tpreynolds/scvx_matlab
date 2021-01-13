function save_data(obj)

dir = '../data/';

scales = obj.get_scales;

filename = strcat(dir,obj.name,'_data');

scvx = saveobj(obj);

save(filename,'scvx','scales');

end

