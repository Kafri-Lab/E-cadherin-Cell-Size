function name_map = fun()
  % Map from keywords in filenames to pretty animal names

  name_map = containers.Map;
  name_map('tiger') = 'Tiger';
  name_map('fruit bat') = 'Fruit Bat';
  name_map('spalax') = 'Blind Mole Rat';
  name_map('wild bore') = 'Wild Bore';
  name_map('horse') = 'Horse';
  name_map('wolf') = 'Grey Wolf';
  name_map('giraf') = 'Giraffe';
  name_map('cow') = 'Cow';
  name_map('cat') = 'Cat';
  name_map('cotton') = 'Cotton Tamarin';
  name_map('kangaroo') = 'Kangaroo';
  name_map('peccary') = 'Peccary';
  name_map('pig') = 'Pig';
  name_map('oryx') = 'Arabian Oryx';
  name_map('nmr') = 'Naked Mole Rat';
  name_map('NMR') = 'Naked Mole Rat';
  name_map('prairie dog') = 'Prairie Dog';
  % name_map('dog') = 'Dog'; % non-unique name (see prairie dog). So we can't use this as is, need to fix
  name_map('zvi') = 'Gazelle';
  name_map('RW') = 'Psammomys';
  name_map('psamon') = 'Psammomys';
  name_map('rat') = 'Black Rat';
  name_map('mouse') = 'Mouse';
  name_map('mou ') = 'Mouse';
  name_map('shrew liv') = 'Shrew';
  name_map('shrew 11') = 'Shrew';
  name_map('human') = 'Human';
  name_map('porcupine') = 'Porcupine';
  name_map('dorban') = 'Porcupine';
  name_map('mon  ') = 'Macaque';
  name_map('grey bat') = 'Grey Bat';
end