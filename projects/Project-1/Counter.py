


class Counter(object):
    
    
    
    
    def __init__(self,Text,kmer_size):
            self.Text = Text
            self.Genome = ['A','C','G','T']
            self.kmer_size = kmer_size
            '''Initiate Frequency Array'''
            self.FrequencyArray = [0 for i in range(4**self.kmer_size)]
            
            
           
    def PatternCount(self,Pattern,verbose=False):
        '''This function counts the given "Pattern" within the supplied "Text":
            Text: DNA string
            Pattern: certain kmer'''
        count = 0
        window_step = 0
        for i in range(len(self.Text) - len(Pattern) + 1):
            '''this counts the sliding window of text'''
            window_step += 1
            
            if verbose:
                print('loop:',i,\
                      'window_step:',window_step,\
                      'count:',count,\
                      "pattern within window:",self.Text[i:i+len(Pattern)])
                
            '''if pattern matches with text within the slidding window then
            count is increased by 1'''  
            if self.Text[i:i+len(Pattern)] == Pattern:
                count += 1 
        return count

    
   
    
    
    def FrequentWords(self , k, verbose=False):
        '''This function finds the frequent kmer of size k by utilizing function "PatternCount()":
            Text: supplied text body of DNA
            k: length of kmer '''
        FreqPatterns = []
        Counts = []
        for i in range(len(self.Text) - k + 1):
            '''Pattern picks up the letter within the sliding window[i:i+k] only'''
            Pattern = self.Text[i:i+k]
            '''Count the pattern frequency within Text'''
            Num = self.PatternCount(self.Text,Pattern)
            '''Add count value in Counts list'''
            Counts.append(Num)
            if verbose:
                print('selected pattern:',Pattern)
                print('pattern num:', Num)
                print('list growth:', Counts)
                print('------------------------')
        '''Find the maximum value of counts'''    
        maxcount = max(Counts)
        for i in range(len(self.Text) - k + 1):
            '''Find out the kmer which has count value equal to maxcount'''
            if Counts[i] == maxcount:
                FreqPatterns.append(self.Text[i:i+k])
                if verbose:
                    print('frequent pattern:',self.Text[i:i+k],\
                          "frequent pattern collection:",FreqPatterns)
        '''Select the unique kmers only'''        
        FreqPatterns = set(FreqPatterns)
        return list(FreqPatterns)
    
    
    
    
    

    def ReverseComplement(self):
        '''initiate a vacant string'''
        rText =""
        '''all letters are replaced by reverse letters'''
        for i in range(len(self.Text)):
            if self.Text[i] == "A":
                  rText = rText+"T"
            elif self.Text[i] == "T":
                  rText  = rText+"A"
            elif self.Text[i] == "G":
                  rText  = rText+"C"
            elif self.Text[i] == "C":
                  rText = rText+"G"
        return rText
           
        
        
        
    def Skew(self, verbose=False):
        '''This function finds the skew score #G-#C of the given Genome'''
        '''SkewList stores current skew for ploting purpose'''
        SkewList = [0]
        '''CurrentSkew keeps changing'''
        CurrentSkew = 0
        '''Iterate through each letter of genome'''
        for i in range(len(self.Text)):
            if self.Text[i] == 'C':
                '''If you meet letter C, CurrentSkew decreases by 1'''
                SkewList.append(SkewList[i] - 1)
                CurrentSkew -= 1
            elif self.Text[i] == 'G':
                '''If you meet letter G, CurrentSkew increases by 1'''    
                SkewList.append(SkewList[i] + 1)
                CurrentSkew += 1
            else:
                '''Do nothing otherwise, just keep appending old CurrentSkew''' 
                SkewList.append(SkewList[i])
            if verbose:
                print(self.Text[i],CurrentSkew)
        return SkewList    
        
        
        
        
        
        
        
    def HammingDistance(self,genome1,genome2):
        '''This function finds the hamming distance by calculating 
        the number of mismatches between two genomes'''
        '''Initiate the score'''
        score = 0
        '''Are genome length equal? send error otherwise'''
        if len(genome1) != len(genome2):
            print("Error")
        else:
            ''' score increases by 1 when mismatch occurs'''
            for i in range(len(genome1)):
                if genome1[i] != genome2[i]:
                    score += 1
        return score     
        
        
        
        
        
        
    def Neighbors(self,Pattern,d, verbose = False):
        '''This function finds the neighbors of a 'Pattern' within hamming distance - 'd' '''
        '''special cases'''
        if d == 0:
            return Pattern
        if len(Pattern) ==1:
            return {'A','C','G','T'}
        '''general case'''
        Neighborhood =  []
        if verbose:
            SuffixNeighbors = self.Neighbors(Pattern[1:],d, verbose=True)
        else:
            SuffixNeighbors = self.Neighbors(Pattern[1:],d)
        if verbose: 
            print("----------------------------")    
            print("SuffixNeighbors:",SuffixNeighbors)
            print("----------------------------") 
        '''Iterate through items in suffixneighbors'''
        for item in SuffixNeighbors:
        
            '''If hmm-distance is below d, add random prefix'''
            if self.HammingDistance(Pattern[1:],item)<d:
                for x in ['A','C','G','T']:
                    Neighborhood.append(x+item)
                if verbose: 
                    print("hamm-dis:",self.HammingDistance(Pattern[1:],item),\
                          "x:",x,\
                          "suffix-item:", item,\
                          "nbr:",x+item)
            else:
                '''otherwise you are not allowed to add random prefix, add original prefix'''  
                Neighborhood.append(Pattern[0]+item)
                if verbose: 
                    print("hamm-dis:",self.HammingDistance(Pattern[1:],item),\
                          "Prefix:",Pattern[0],\
                          "suffix-item:", item,\
                          "nbr:",Pattern[0]+item )
        return Neighborhood  
    
    
    
    def ApproximatePatternCount(self, kmer, hamming_score,verbose=False):
        '''This function finds frequency of a pattern with
        neighbours within 'hamming_score' '''
        '''Initiate Approximate Count'''
        Count = 0
        '''Set the slidding window'''
        for i in range(len(self.Text) - len(kmer) + 1):
            '''kmer within the window is Pattern'''
            Pattern = self.Text[i:i+len(kmer)]
            Score = self.HammingDistance(kmer,Pattern)
            '''If hmm-dist is below hamming_score, keep increasing  Approximate Count'''
            if Score <= hamming_score:
                Count += 1
                if verbose:
                    print("window-run::",i,"to",i+len(kmer),\
                          "Pattern:", genome[i:i+len(kmer)],\
                          "count:",Count)
        return Count
    
    
    def SymbolToNumber(self, symbol):
        for s in enumerate(self.Genome):
            if s[1]==symbol:
                return s[0]
    
    
    
    
    def NumberToSymbol(self,n):
        if n==0:
            return "A"
        elif n==1:
            return "C"
        elif n==2:
            return "G"
        else:
            return "T"
        
        
        
    def PatternToNumber(self,Pattern,verbose=False):
    
        '''This function returns a index of a 'Pattern' using recursive algorithm'''
            
        if len(Pattern)==0:
        
            return 0
    
        else:
    
            '''last letter is a symbol'''
            Symbol = Pattern[-1]
    
            '''all but last are prefix pattern'''
            Prefix = Pattern[0:-1]
    
    
            if verbose:
                print("prefixPattern:",Prefix,\
                  ", Symbol:",Symbol,\
                  ", prefixIndex:",self.PatternToNumber(Prefix),\
                  ", SymbolIndex:",self.SymbolToNumber(Symbol),\
                  "=>",\
                  "4*(",self.PatternToNumber(Prefix),\
                  ")+(",self.SymbolToNumber(Symbol),\
                  ") = ",4*self.PatternToNumber(self.Prefix)+self.SymbolToNumber(Symbol)  )
            
        
    
            if verbose:
                '''to keep printing in every cycle''' 
                return 4*self.PatternToNumber(Prefix,verbose=True)+self.SymbolToNumber(Symbol)
            else:
                '''verbose=False id defult here'''
                return 4*self.PatternToNumber(Prefix)+self.SymbolToNumber(Symbol)  
        
        
    def NumberToPattern(self,Index,PatternLength,verbose=False):
    
        '''This function returns Pattern for given Index'''
        
        if PatternLength ==1:
            return self.NumberToSymbol(Index)
        
        else:
            
            '''Index = 4*quotient+ reminder'''
            prefixIndex = int(Index/4)
            Reminder = Index%4
            
            
            '''Convert reminder to Symbol'''
            Symbol = self.NumberToSymbol(Reminder)
            
           
            if verbose:
                print("prefixIndex:",prefixIndex,\
                      ", Reminder:",Reminder,
                      ", Symbol:",Symbol,\
                      "=>",Index,"=",\
                      "4*",prefixIndex,\
                      "+",Reminder)
                
                
            if verbose:
                return self.NumberToPattern(prefixIndex,PatternLength-1, verbose =True) + Symbol
            else:
                return self.NumberToPattern(prefixIndex,PatternLength-1) + Symbol
            
            
            
    def computingFrequrncies(self,verbose=False):
    
      
        '''Iterate through sliding window'''
        for i in range(len(self.Text)-self.kmer_size):
        
        
            '''Pattern is the text within the sliding window'''
            Pattern = self.Text[i:i+self.kmer_size]
        
        
            '''Find the index of the selected pattern'''
            Index = self.PatternToNumber(Pattern)
        
            '''Increase frequency of the pattern in FrequencyArray'''
            self.FrequencyArray[Index] = self.FrequencyArray[Index]+1
        
        
        
            if verbose:
                print("Pattern:", Pattern,\
                      "Index:",Index, \
                      "Increasing FrequencyArray:", self.FrequencyArray[Index])
            
            
        return FrequencyArray
    