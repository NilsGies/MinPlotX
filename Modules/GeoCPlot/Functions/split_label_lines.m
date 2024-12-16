function string_skalar = split_label_lines(label,max_label_size)

label_text=label;
if iscell(label)
label_size=numel(label{:});

    string_skalar=label;

else
    label_size=numel(label);

    string_skalar={label};
end
if label_size > max_label_size
    string_skalar = {}; % cell array of strings that will contain the final xlabel String
    splitted_label_text = split(label_text); % split the text into words
    n_words = numel(splitted_label_text); % count the number of words

    label_size = 0; % xlabel width in pixels
    n_words_new = 0; % number of words on the current line of the label
    word_idx = 1; % index of the next word to be added to the label
    last_word_size = 0; % number of characters in the last word added to the label
    while word_idx <= n_words
        if n_words_new == 0
            current_line = splitted_label_text{word_idx};
            last_word_size = numel(splitted_label_text{word_idx});
        else
            current_line = [current_line ' ' splitted_label_text{word_idx}];
            last_word_size = numel(splitted_label_text{word_idx})+1;
        end
        n_words_new = n_words_new+1;
        word_idx = word_idx+1;


        % and get the xlabel's width in pixels
        label_size=numel(current_line);

        if label_size > max_label_size
            if n_words_new == 1
                string_skalar{end+1} = current_line;
                current_line = '';
                n_words_new = 0;
            else
                string_skalar{end+1} = current_line(1:end-last_word_size);
                current_line = current_line(end-last_word_size+2:end);
                n_words_new = 1;
            end
        end
    end
    if ~isempty(current_line)
        string_skalar{end+1} = current_line;
    end
end
