import pandas as pd

def RelativeAbundance(OTU_table_full_taxonomy, Query_taxonomy_level):
    """
    Reurn relative abundance dataframe at different taxonomic level.
    OTU_table_full_taxonomy: dataframe, columns are OTU ID, Sample name, taxonomy
    Query_taxonomy_level: str, one of ['k', 'p', 'c', 'o', 'f', 'g', 's']
    """
    #复制otu_table dataframe，并把OTU ID设置为index
    otu_table_abundance = OTU_table_full_taxonomy.copy(deep=True)
    otu_table_abundance.set_index("#OTU ID", inplace=True)

    #拆分taxonomy数据为字典
    def parse_taxonomy(taxonomy_str):
        taxonomy_levels = ['k', 'p', 'c', 'o', 'f', 'g', 's'] 
        taxonomy_dict = {level: 'unclassified' for level in taxonomy_levels}
        items = taxonomy_str.split(',')
        for item in items:
            key, value = item.split(':')
            if key.strip() in taxonomy_dict:
                taxonomy_dict[key.strip()] = value.strip()
        return taxonomy_dict
    
    taxonomy_df = otu_table_abundance['taxonomy'].apply(parse_taxonomy).apply(pd.Series)
    otu_table_abundance = otu_table_abundance.drop(columns='taxonomy').join(taxonomy_df)

    #不同分类水平的otu_table
    taxonomic_otu_table = otu_table_abundance.groupby(Query_taxonomy_level).sum()

    #计算taxonomic_otu_table相对丰度
    #columns: samples
    #rows: OTU relative abundance
    taxonomic_relative_abundance = taxonomic_otu_table.div(taxonomic_otu_table.sum(axis=0), axis=1)

    #转置taxonomic otu table
    #columns: OTu
    # rows: samples
    taxonomic_relative_abundance = taxonomic_relative_abundance.T
    return taxonomic_relative_abundance

if __name__ == '__main__':
    main()