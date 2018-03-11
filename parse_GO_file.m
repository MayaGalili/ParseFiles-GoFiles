%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script goal is to take GO file and turn it into anotation
% structure. The structure is needed for SAFE analysis (Systematic
% functional annotation and visualization) and other large scale genetic
% analysis.
% 
% INPUT:	GO_file_nm: file path
%
% OUTPUT: 	GO_mat: GXT size binary matrix. [GO_mat(g,t)=1 if gene g known as connected to term t]
%			gene_set: G size vector with genes names. sorted alphabeticly.
%			terms_set: T size vector with terms names. sorted alphabeticly.
%
% @ Maya Galili. Dec 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [go] = parse_GO_file(GB_file_nm)

% init vars
go =[]; % create the GO structure
go.modified_on=date;
go.description =[ 'GO structure contains: '...
    '{term_ids: GO terms list} ' ...
    '{term_names: GO terms description} ' ...
    '{orfs: ORFs annotated to terms} ' ...
    '{term2orf: binary matrix of genes mapped to ORFs}'];
	
% load and read the file
goGenes = goannotread(GB_file_nm);
gene_list ={goGenes.TaxonomyID}';
terms_list = {goGenes.GOterm}';
terms_list=regexprep(terms_list, 'GO:', '');
terms_list=cellfun(@(x) str2double(x), terms_list,'UniformOutput' , false);
terms_list = cell2mat(terms_list);

% find Term list and sort by names
unique_go_list = sort(unique(terms_list));
unique_go_list(isnan(unique_go_list)) = [];
[unique_go_list, term_names_list] = from_term_ids_to_term_names(unique_go_list);

% find Gene list and sort by names
unique_geneORF_list = sort(unique(gene_list));

% create the GO matrix
termOrf_mat = get_term2orf_mat(unique_geneORF_list, unique_go_list,...
    gene_list, terms_list);

% remove terms with no genes
terms_to_rmv_idx=find(sum(termOrf_mat,2)==0);
termOrf_mat(terms_to_rmv_idx,:) = [];
unique_go_list(terms_to_rmv_idx)=[];

go.orfs = unique_geneORF_list;  % B string cell array
go.term2orf = termOrf_mat;  % M=AXB binary (uint8)
go.term_ids=unique_go_list;   % A double array
go.term_names=term_names_list; % A string cell array

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fonction goal: from GO ID to Term discription
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [terms_ids, terms_names] = from_term_ids_to_term_names(terms_ids)

% load all GO terms
GO_anot = geneont('LIVE', true); 

% find listed terms dicription
subontology = GO_anot(terms_ids);
terms_names = get(subontology.Terms,{'name'}); % A string cell array
terms_ids = get(subontology.Terms,{'id'}); % A double array
terms_ids = cell2mat(terms_ids);
end
