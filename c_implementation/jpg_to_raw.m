clc ; clear ll; close all;

jpg_folder = ('D:\vot2014\ball\color\*.jpg');
folder = ('D:\vot2014\ball\color\');
file_list=dir(jpg_folder);
file_num=size(file_list,1);

for i=1:file_num
  list1 = {file_list.name};  
  list2 = list1';
 
  file_name = strcat(folder,list2{i});
  img=imread(file_name);%读取图像
  img_1=img(:);
  name_s = strsplit(file_name,'.');
  
  txt_name = strcat(name_s{1},'.txt');
  fid=fopen(txt_name,'w');%存为txt
  fprintf(fid, '%d\n',img_1);%注意将img转置
  fclose(fid);
  
end  
endl = 1 

    

