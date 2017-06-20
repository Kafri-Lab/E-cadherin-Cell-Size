function animal_name = img_name_to_animal_name(img_name,name_map)
  key = keys(name_map); % get keys
  animal_name = 'not found';

  for n=1:name_map.Count
    if strfind(img_name,key{n})
      animal_name = name_map(key{n});
      break
    end
  end

end
