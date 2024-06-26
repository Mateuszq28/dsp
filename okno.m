function okno = okno(typ, rozmiar);
  
  switch (typ)
    case (1)
      printf('Mialo byc kazde z wyjatkiem tego!\n')
    case (2)
      okno = bartlett(rozmiar);
    case (3)
      okno = blackman(rozmiar);
    case (4)
      okno = hamming(rozmiar);
    case (5)
      okno = hanning(rozmiar);
    otherwise
      printf("Niewlasciwy typ okna \n")
  end
  
  okno = round(okno*32767);                      %32768 to za du¿o
  
  fileID = fopen('Jakas_nazwa_tablicy_okna.h','w');
  fprintf(fileID, 'const int Jakas_nazwa_tablicy_okna[%d] = {\n', rozmiar);
  for i = 1:(rozmiar)
    if numel(num2str(okno(i))) == 4
      if rem(i, 16) == 0
        fprintf(fileID, ' %d, \n', okno(i));
      else    
        fprintf(fileID, ' %d, ', okno(i));
      end
    else
      if rem(i, 16) == 0
        fprintf(fileID, '%d,\n', okno(i));
      else    
        fprintf(fileID, '%d, ', okno(i));
      end
    end
  end
    
  fprintf(fileID,'};\n', okno(i));
  fclose(fileID);
  
end