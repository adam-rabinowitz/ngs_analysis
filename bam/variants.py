import os
import pysam
import collections
import multiprocessing

class basecall(object):
    
    def __init__(
            self, bam
        ):
        self.bam = bam
        bam = pysam.AlignmentFile(bam)
        self.names = collections.OrderedDict(
            {bam.gettid(x):x for x in bam.references})
        self.lengths = {x:y for x, y in zip(bam.references, bam.lengths)}
        bam.close()
    
    def extract_reads(
            self, pipe_out, intervals, remove_dup = False, remove_sec = False,
            map_quality = 0
        ):
        ''' Function to read FASTQ files and send data down pipe '''
        print('extract_reads: {}'.format(os.getpid()))
        # Set flag
        flag = 4
        if remove_sec:
            flag += 256
        if remove_dup:
            flag += 1024
        # Loop through reads in chromsome
        bam = pysam.AlignmentFile(self.bam)
        for chrom, start, end in intervals:
            pipe_out.send((chrom, start, end))
            for read in bam.fetch(chrom, start, end):
                if read.flag & flag:
                    continue
                if read.mapping_quality < map_quality:
                    continue
                read_data = {
                    'sequence' : list(read.query_alignment_sequence),
                    'basequal' : list(read.query_alignment_qualities),
                    'cigar' : read.cigar,
                    'mapqual' : read.mapping_quality,
                    'forward' : not(read.is_reverse),
                    'refstart' : read.reference_start,
                    'refend' : read.reference_end
                }
                pipe_out.send(read_data)
        pipe_out.close()
        bam.close()
    
    def extract_reads_process(
            self, intervals, remove_dup = False, remove_sec = False,
            map_quality = 0
        ):
        connRecv, connSend = multiprocessing.Pipe(False)
        process = multiprocessing.Process(
            target = self.extract_reads,
            args = (connSend, intervals, remove_dup, remove_sec, map_quality)
        )
        process.start()
        connSend.close()
        return((connRecv, process))

    def extract_bases(
            self, pipe_in, pipe_out, ignore_indel = False,
            group_del = False, base_quality = 0
        ):
        print('extract_bases: {}'.format(os.getpid()))
        while True:
            try:
                read_data = pipe_in.recv()
            except EOFError:
                break
            # Unpack data
            try:
                cigar = read_data.pop('cigar')
                sequence = read_data.pop('sequence')
                basequal = read_data.pop('basequal')
                refend = read_data.pop('refend')
            except AttributeError:
                pipe_out.send(read_data)
                ref_chrom, ref_start, ref_end = read_data
                continue
            # Remove clipping and terminal indels from cigar
            cigar = [x for x in cigar if x[0] not in (4,5)]
            while cigar:
                if cigar[0][0] in (1,2):
                    sequence = sequence[cigar[0][1]:]
                    basequal = basequal[cigar[0][1]:]
                    cigar = cigar[1:]
                else:
                    break
            while cigar:
                if cigar[-1][0] in (1,2):
                    sequence = sequence[:-cigar[-1][1]]
                    basequal = basequal[:-cigar[-1][1]]
                    cigar = cigar[:-1]
                else:
                    break
            if not cigar:
                continue
            # Extract quality
            base_list = []
            read_index = 0
            ref_index = read_data['refstart']
            for operation, operation_len in cigar:
                # Process match operations
                if operation == 0:
                    for _ in range(operation_len):
                        if basequal[read_index] > base_quality:
                            base_list.append((
                                ref_index,
                                sequence[read_index],
                                basequal[read_index],
                            ))
                        read_index += 1
                        ref_index += 1
                # Process insertion operations
                elif operation == 1:
                    if ignore_indel:
                        read_index += operation_len
                    else:
                        # Adjust indices to include previous base
                        read_index -= 1
                        ref_index -= 1
                        operation_len += 1
                        # Create insertion slice and calculate mean quality
                        insert_slice = slice(read_index,
                            read_index + operation_len)
                        mean_quality = (sum(basequal[insert_slice]) /
                            operation_len)
                        # Replace previous base call with insertion
                        if mean_quality >= base_quality:
                            insert_sequence = ''.join(sequence[insert_slice])
                            base_list.pop()
                            base_list.append(
                                (ref_index, insert_sequence, mean_quality)
                            )
                        # Adjust read and reference indices
                        read_index += operation_len
                        ref_index += 1
                # Process deletion operations
                elif operation == 2:
                    if ignore_indel:
                        ref_index += operation_len
                    else:
                        # Calculate mean quality of bases surrounding deletion
                        mean_quality = sum(basequal[read_index - 1 : read_index + 1]) / 2
                        if mean_quality >= base_quality:
                            # Add base call for grouped delet
                            if group_del:
                                delete_sequence = '-' * operation_len
                                base_list.append(
                                    (ref_index, delete_sequence, mean_quality)
                                )
                                ref_index += operation_len
                            # Add base calls for non grouped deletion
                            else:
                                for _ in range(operation_len):
                                    base_list.append(
                                        (ref_index, '-', mean_quality)
                                    )
                                    ref_index += 1
                # Raise error if other qualoties found
                else:
                    raise ValueError('Unexpected cigar {}'.format(cigar))
            # Check processing
            if (refend != ref_index
                or read_index != len(sequence)):
                raise ValueError('Error processing read')
            # Add data to read data dictionary
            base_list = [x for x in base_list if ref_start <= x[0] <= ref_end]
            if base_quality:
                base_list = [x for x in base_list if x[2] >= base_quality]
            read_data['baselist'] = base_list
            pipe_out.send(read_data)
        # Clean up
        pipe_in.close()
        pipe_out.close()

    def extract_bases_process(
            self, intervals, remove_dup = False, remove_sec = False,
            map_quality = 0, base_quality = 0, ignore_indel = False,
            group_del = True
        ):
        # Create process to read pipe
        read_pipe, read_process = self.extract_reads_process(
            intervals = intervals, remove_dup = remove_dup,
            remove_sec = remove_sec, map_quality = map_quality
        )
        # Create process to extract reads
        conn_recv, conn_send = multiprocessing.Pipe(False)
        process = multiprocessing.Process(
            target = self.extract_bases,
            args = (read_pipe, conn_send, ignore_indel, group_del,
                base_quality)
        )
        process.start()
        conn_send.close()
        # Return data
        return(conn_recv, [read_process, process])


        base_pipe, base_process = self.extract_bases(
            self, pipe_in, pipe_out, ignore_indel = False,
            group_del = False, base_quality = 0
        )
        # Start read extraction process


    def base_calls_meanmap(
            self, pipe_in, pipe_out, map_quality = 0, base_quality = 0
        ):
        print('base_calls: {}'.format(os.getpid()))
        # Create variables for storage
        while True:
            # Extract data from pipe or flush current data
            try:
                read_data = pipe_in.recv()
            except EOFError:
                for position in base_dict:
                    base_data = base_dict.pop(position)
                    mean_map = base_data[0][0] / base_data[0][1]
                    pipe_out.send(
                        (chrom, position, mean_map, base_data[1])
                    )
                break
            # Extract new interval or reference start from read
            try:
                new_start = read_data['baselist'][0][0]
            except TypeError:
                base_dict = collections.OrderedDict()
                chrom = read_data[0]
                current_start = 0
                continue
            # Extract positions left of reference start
            if new_start < current_start:
                raise ValueError('Is BAM sorted?')
            for position in base_dict:
                if position < new_start:
                    base_data = base_dict.pop(position)
                    mean_map = base_data[0][0] / base_data[0][1]
                    pipe_out.send(
                        (chrom, position, mean_map, base_data[1])
                    )
                else:
                    break
            # Add data to new start
            for position, call, basequal in read_data['baselist']:
                if position not in base_dict:
                    base_dict[position] = [[0, 0], collections.Counter()]
                # Append
                base_dict[position][0][0] += read_data['mapqual']
                base_dict[position][0][1] += 1
                if (basequal >= base_quality
                    and read_data['mapqual'] >= map_quality):
                    base_dict[position][1][call] += 1
            # Reset start
            current_start = new_start
        # Clean up
        pipe_in.close()
        pipe_out.close()
    
    def base_calls(
            self, pipe_in, pipe_out, map_quality = 0, base_quality = 0
        ):
        print('base_calls: {}'.format(os.getpid()))
        # Create ordered dictionary for storage
        base_dict = collections.OrderedDict()
        # Collate read data
        while True:
            # Extract data from pipe or flush current data
            try:
                read_data = pipe_in.recv()
            except EOFError:
                break
            # Extract new interval or reference start from read
            try:
                new_start = read_data['refstart']
            except TypeError:
                # Empty base dictionary and initialise variables
                for position in base_dict:
                    base_data = base_dict.pop(position)
                    pipe_out.send((chrom, position, base_data))
                chrom = read_data[0]
                current_start = 0
                continue
            # Extract positions left of reference start
            if new_start < current_start:
                raise ValueError('Is BAM sorted?')
            for position in base_dict:
                if position < new_start:
                    base_data = base_dict.pop(position)
                    pipe_out.send((chrom, position, base_data))
                else:
                    break
            # Skip poorly mapped reads
            if read_data['mapqual'] < map_quality:
                continue
            # Add base calls to dictionary
            for position, call, basequal in read_data['baselist']:
                if basequal >= base_quality:
                    try:
                        base_dict[position][call] += 1
                    except KeyError:
                        base_dict[position] = collections.defaultdict(int)
                        base_dict[position][call] += 1
            # Reset start
            current_start = new_start
        # Clean up
        pipe_in.close()
        for position in base_dict:
            base_data = base_dict.pop(position)
            pipe_out.send((chrom, position, base_data))
        pipe_out.close()

    def base_call_process(
            self, intervals, remove_dup = False, remove_sec = False,
            ignore_indels = False, group_del = True, map_quality = 0,
            base_quality = 0, mean_map = False
        ):
        conn_recv, conn_send = multiprocessing.Pipe(False)
        # Create process to extract reads
        if mean_map:
            base_pipe, process_list = self.extract_bases_process(
                intervals = intervals, remove_dup = remove_dup,
                remove_sec = remove_sec, map_quality = 0, base_quality = 0)
            call_process = multiprocessing.Process(
                target = self.base_calls_meanmap,
                args = (base_pipe, conn_send, map_quality, base_quality)
            )
        else:
            base_pipe, process_list = self.extract_bases_process(
                intervals = intervals, remove_dup = remove_dup,
                remove_sec = remove_sec, map_quality = map_quality,
                base_quality = base_quality)
            call_process = multiprocessing.Process(
                target = self.base_calls,
                args = (base_pipe, conn_send, map_quality, base_quality)
            )
        # Start process
        call_process.start()
        process_list.append(call_process)
        conn_send.close()
        return(conn_recv, process_list)

x = basecall('/farm/scratch/rs-bio-lif/rabino01/Elza/bamFiles/310123-NORM.chr19.bam')
pipe, process = x.base_call_process([('19', 30716400, 30716700)], mean_map = False)
counter = 0
while True:
    try:
        data = pipe.recv()
        counter += 1
        print data
    except EOFError:
        break
print counter
