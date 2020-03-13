function [approx , err_abs] = bissec(f , x0 , x1 , nb_it_max , tol_rel , file_name)
% BISSEC	M?thode de la bissection pour la r?solution f(x) = 0
%			pour f: R -> R
%
% Syntaxe: [approx , err_abs] = bissec(f , x0 , x1 , nb_it_max , tol_rel , file_name)
%
% Arguments d'entr?e
%	f			-	String ou fonction handle sp?cifiant la fonction
%					non-lin?aire
%	x0			-	1?re approximation initiale 
%	x1			-	2?me approximation initiale 
%	nb_it_max	-	Nombre maximum d'it?rations 
%	tol_rel		-	Tol?rance sur l'approximation de l'erreur relative
%	file_name	-	(Optionnel) Nom du fichier (avec l'extension .txt) dans
%					lequel sera	?crit les r?sultats de l'algorithme
%
% Arguments de sortie
%	approx		-	Vecteur colonne de taille nb_iter contenant les 
%					it?rations
%	err_abs		-	Vecteur colonne de dimension nb_iter contenant les
%					erreurs absolues
%
% Exemples d'appel
%	[ approx , err_abs ] = bissec( 'my_fct_nl' , 3 , 4 , 100 , 1e-9 , 'resul_bissec.txt')
%	[ approx , err_abs ] = bissec( @(x) x.^2-10 , 3 , 4 , 100 , 1e-9 , 'resul_bissec.txt')




%%  V?rification de la fonction
if isa(f,'char')
	fct		=	str2func(f);
	is_fct_file = true;
elseif isa(f,'function_handle')
	fct		=	f;
	is_fct_file = false;
else
	error('L''argument f n''est pas un string ni un function_handle')
end

%% V?rification du nb de composantes de l'approximation initiale et de f
if ~isnumeric(x0) || ~isscalar(x0)
	error('L''approximation initiale x0 n''est pas un scalaire')
elseif ~isnumeric(x1) || ~isscalar(x1)
	error('L''approximation initiale x1 n''est pas un scalaire')
end

try 
	fct(x0);
catch ME
	if strcmp(ME.identifier,'MATLAB:UndefinedFunction')
		error('La fonction f n''est pas dans le r?pertoire courant')
	elseif strcmp(ME.identifier,'MATLAB:badsubscript')
		error('Le fonction f ne retourne pas un scalaire')
	elseif strcmp(ME.identifier,'MATLAB:minrhs')
		error('La fonction f doit prendre seulement 1 argument en entr?e')
	else 
		rethrow(ME)
	end
end

if ~isnumeric(fct(x0)) || ~isscalar(fct(x0)) 
	error('Le fonction f ne retourne pas un scalaire')
elseif fct(x0)*fct(x1)>0
	warning('La condition f(x0)*f(x1)<0 n''est pas respect?e.\nArr?t de l''algorithme%\n',[])
	approx	=	[x0;x1];
	err_abs =	inf(2,1);
	return
elseif fct(x0)==0
	approx	=	x0;
	err_abs	=	0;
	return
elseif fct(x1)==0
	approx	=	x1;
	err_abs	=	0;
	return
end

%% V?rification du fichier output
if nargin == 6 && ~isa(file_name,'char')
	error('Le nom du fichier des r?sultats doit ?tre de type string')
end

%% Initialisation des matrices app et err
app			=	nan(nb_it_max,1);
err_rel		=	inf(nb_it_max,1);
arret		=	false;


%% M?thode de la bissection
for t=1:nb_it_max-1
	
	if t==1
		x_gauche	=	min([x0,x1]);
		x_droite	=	max([x0,x1]);
	else
		if f_gauche*f_milieu < 0
			x_droite	=	x_milieu;
		elseif f_droite*f_milieu < 0
			x_gauche	=	x_milieu;
		else
			warning('Probl?me avec la fonction f.\nArr?t de l''algorithme%\n',[])
			break
		end
	end
	
	x_milieu	=	(x_gauche + x_droite)/2;
	app(t)		=	x_milieu;
	
	if t==1
		if fct(app(t)) == 0
			arret	=	true;
			break
		end
	else
		err_rel(t-1)	=	abs(app(t)-app(t-1))/(abs(app(t)) + eps);
		
		if (err_rel(t-1) <= tol_rel) || (fct(app(t)) == 0)
			arret	=	true;
			break
		end
	end

	f_gauche	=	fct(x_gauche);
	f_droite	=	fct(x_droite);
	f_milieu	=	fct(x_milieu);
	
end

nb_it	=	t;
approx	=	app(1:nb_it);
err_abs	=	inf(nb_it,1);

if arret
	err_abs		=	abs(approx(end) - approx);
else
	warning('La m?thode de la bissection n''a pas converg?e')
end

% ?criture des r?sultats si fichier pass? en argument
if nargin == 6
	output_results(file_name , fct , is_fct_file , nb_it_max , ...
					tol_rel , x0 , x1 , approx , err_abs , arret)
end

end


function [] = output_results(file_name , fct , is_fct_file , it_max , ...
							tol_rel , x0 , x1 , x , err , status)
						 
	fid		=	fopen(file_name,'w');
	fprintf(fid,'Algorithme de la bissection\n\n');
	fprintf(fid,'Fonction dont on cherche les racines:\n');
	if is_fct_file
		fprintf(fid,'%s\n\n',fileread([func2str(fct),'.m']));
	else
		fprintf(fid,'%s\n\n',func2str(fct));
	end
	fprintf(fid,'Arguments d''entr?e:\n');
	fprintf(fid,'    - Nombre maximum d''it?rations: %d\n',it_max);
	fprintf(fid,'    - Tol?rance relative: %6.5e\n',tol_rel);
	fprintf(fid,'    - Approximation initiale x0: %16.15e\n',x0);
	fprintf(fid,'    - Approximation initiale x1: %16.15e\n\n',x1);
	
	if status
		fprintf(fid,'\nStatut: L''algorithme de la bissection a converg? en %d it?rations\n\n',length(x));
	else
		fprintf(fid,'\nStatut: L''algorithme de la bissection n''a pas converg?\n\n');
	end
	fprintf(fid,'#It            x              Erreur absolue\n');
	fprintf(fid,'%3d   %16.15e   %6.5e\n',[reshape(1:length(x),1,[]);reshape(x,1,[]);reshape(err,1,[])]);
	fclose(fid);

end